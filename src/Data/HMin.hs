{-# LANGUAGE FlexibleContexts, TypeOperators, Rank2Types, TypeFamilies #-}

module Data.HMin where

import Data.MemoTrie
import Data.AdditiveGroup
import Data.VectorSpace
import Data.Maclaurin
import Data.Basis
import Data.LinearMap
import Control.Arrow

-- convenience type and eq operator
type Equation a b = (a -> b) -> b -> a -> b

infixl 5 ~=
(~=) :: Num b => (a -> b) -> b -> a -> b
f ~= c = \x -> f x - c


eqAnd :: Num c => (a -> c) -> (b -> c) -> (a, b) -> c
eqAnd f1 f2 (x, y) = sqr (f1 x) + sqr (f2 y)

eqAnd' :: Num c => (a -> c) -> (a -> c) -> a -> c
eqAnd' f1 f2 x = (f1 `eqAnd` f2) (dup x)
                    where dup v = (v, v)

sqr :: Num b => b -> b
sqr x = x*x

sqrD :: (Num b, VectorSpace b, b ~ Scalar b, b ~ (Scalar a), HasBasis a, HasTrie (Basis a)) => (a :> b) -> (a :> b)
sqrD = sqr


solve :: (Show a, Num b, Ord b, AdditiveGroup a) => (a :~> b) -> ((a :> b) -> Bool) -> ((a :> b) -> a) -> a -> a
solve eq stopfunc stepfunc start = let delta = eq start in
                                if stopfunc delta
                                    then start
                                    else solve eq stopfunc stepfunc (start ^+^ stepfunc delta)

-- TODO
-- need this
-- https://en.wikipedia.org/wiki/Pushforward_(differential)
gradStepFunc :: (HasBasis a, HasTrie (Basis a), VectorSpace b, Num b, Scalar a ~ b, Scalar b ~ b) => a -> (a :> b) -> a
gradStepFunc gammas delta = negateV $ dV gammas (sqrD delta)

-- TODO double check this please
-- the derivative in the V direction with "f," a function to get us
-- from a :> b to a...
dV :: (Scalar a ~ b, Scalar b ~ b, HasBasis a, HasTrie (Basis a), VectorSpace b) => a -> (a :> b) -> a
dV dv dfx = recompose . map (\(v, s) -> (v, (^* s) . powVal $ derivAtBasis dfx v)) $ decompose dv


equ' :: ((Double, Double) :~> Double)
equ' = sqr fstD + sqr sndD ~= 4

equ'' :: ((Double, Double) :~> Double)
equ'' = fstD + 2*sndD ~= 2

-- this seems to requre the list of basis values...
{-
instance (VectorSpace b, HasBasis a, HasTrie (Basis a), Scalar a ~ Scalar b) => HasBasis (a -> b) where
    type Basis (a -> b) = (Basis a) -> b

    -- TODO
    -- do we really need Scalar a ~ Scalar b?
    -- lapply requires it, so seems possible.
    -- move to linear maps...
    basisValue f = sumV . map (\(b, s) -> s *^ f b) . decompose

    -- TODO
    -- I think I need the list of basis vectors (possibly infinite) here...
    -- Data.VectorSpace.project???
    decompose f = map (first (f . basisValue))
-}

-- https://stackoverflow.com/questions/9313994/derivative-towers-and-how-to-use-the-vector-space-package-haskell

diff :: (Double :~> (Double,Double,Double) ) -> (Double :~> (Double,Double,Double))
diff g = \x -> (atBasis (derivative (g x)) ())

diff' :: Either () () -> ((Double, Double) :~> Double) -> ((Double, Double) :~> Double)
diff' b g = \(x,y) -> derivAtBasis (g (x,y)) b

eval :: (a :~> b) -> a -> b
eval g x = powVal (g x)

f :: Double :~> (Double,Double,Double)
f x = tripleD (pureD 0, pureD 1, (2*idD) x)

f' :: (Double, Double) :~> Double
f' xy = fstD xy + (sndD*2) xy
