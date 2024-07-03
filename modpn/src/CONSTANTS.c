/* This file is part of the MODPN library

    MODPN is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    MODPN is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Foobar.  If not, see <http://www.gnu.org/licenses/>.

    Copyright, Marc Moreno Maza <moreno@csd.uwo.ca>
*/


/* Authors: Xin Li <xli96@csd.uwo.ca>, Marc Moreno Maza <moreno@csd.uwo.ca> */
/* Copyright (c) 2009 by Marc Moreno Maza.  All rights reserved             */
#include "CONSTANTS.h"

/* Variable Table */
char letters[27] = {'0','a','b','c','d','e','f','g','h','i','j','k',
    'l','m','n', 'o','p','q','r','s','t','u','v','w','x','y','z'};

/* The size of machine integer */
#define bitmod 32
sfixn BASE = bitmod;
sfixn BASEHALF = bitmod / 2;
sfixn BASE_1 = bitmod - 1;

/* Fourier Prime Table */
sfixn FP[100] = {
    962592769, 957349889, 950009857, 943718401, 940572673, 938475521,
    935329793, 925892609, 924844033, 919601153, 918552577, 913309697,
    907018241, 899678209, 897581057, 883949569, 880803841, 862978049,
    850395137, 833617921, 824180737, 818937857, 802160641, 800063489,
    799014913, 786432001, 770703361, 754974721, 745537537, 740294657,
    718274561, 715128833, 710934529, 683671553, 666894337, 655360001,
    648019969, 645922817, 639631361, 635437057, 605028353, 597688321,
    595591169, 581959681, 576716801, 531628033, 493879297, 469762049,
    468713473, 463470593, 459276289, 447741953, 415236097, 409993217,
    399507457, 387973121, 383778817, 377487361, 361758721, 359661569,
    347078657, 330301441, 311427073, 305135617, 290455553, 274726913,
    270532609, 257949697, 249561089, 246415361, 230686721, 221249537,
    211812353, 204472321, 199229441, 186646529, 185597953, 169869313,
    167772161, 163577857, 158334977, 155189249, 147849217, 141557761,
    138412033, 136314881, 132120577, 120586241, 113246209, 111149057,
    104857601, 101711873, 81788929,  70254593,  69206017,  28311553,
    26214401,  23068673,  13631489,  7340033};

/* global flag to set if cuda is enabled or not */
sfixn CUDA_TAG = 1;

