OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.090156) q[0];
sx q[0];
rz(-0.12026726) q[0];
sx q[0];
rz(-0.076844849) q[0];
rz(-0.85637561) q[1];
sx q[1];
rz(-2.0937347) q[1];
sx q[1];
rz(0.82077789) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4030006) q[0];
sx q[0];
rz(-1.101491) q[0];
sx q[0];
rz(-1.4856536) q[0];
x q[1];
rz(-1.7834383) q[2];
sx q[2];
rz(-0.66591381) q[2];
sx q[2];
rz(0.41541157) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.1917853) q[1];
sx q[1];
rz(-2.2151183) q[1];
sx q[1];
rz(-0.61442805) q[1];
x q[2];
rz(0.28750136) q[3];
sx q[3];
rz(-0.56351065) q[3];
sx q[3];
rz(-2.7873791) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.0277752) q[2];
sx q[2];
rz(-0.57089388) q[2];
sx q[2];
rz(1.6119733) q[2];
rz(-2.2714603) q[3];
sx q[3];
rz(-1.2711997) q[3];
sx q[3];
rz(2.1421656) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40365264) q[0];
sx q[0];
rz(-1.0920748) q[0];
sx q[0];
rz(-1.8988443) q[0];
rz(-0.87951648) q[1];
sx q[1];
rz(-2.1915235) q[1];
sx q[1];
rz(1.3290728) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91956988) q[0];
sx q[0];
rz(-2.0139512) q[0];
sx q[0];
rz(-0.65624563) q[0];
rz(-pi) q[1];
x q[1];
rz(0.7793541) q[2];
sx q[2];
rz(-2.4930304) q[2];
sx q[2];
rz(0.88403242) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.14790598) q[1];
sx q[1];
rz(-0.68274409) q[1];
sx q[1];
rz(-2.3629689) q[1];
rz(-pi) q[2];
rz(-0.69345052) q[3];
sx q[3];
rz(-1.4458522) q[3];
sx q[3];
rz(3.0378064) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.8326412) q[2];
sx q[2];
rz(-1.2244886) q[2];
sx q[2];
rz(1.857117) q[2];
rz(-2.4744611) q[3];
sx q[3];
rz(-2.7693373) q[3];
sx q[3];
rz(-2.4767806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2892147) q[0];
sx q[0];
rz(-2.7279655) q[0];
sx q[0];
rz(-1.941823) q[0];
rz(1.2780227) q[1];
sx q[1];
rz(-2.8646902) q[1];
sx q[1];
rz(-0.75495458) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2374868) q[0];
sx q[0];
rz(-1.5747935) q[0];
sx q[0];
rz(-2.0275751) q[0];
rz(0.16904449) q[2];
sx q[2];
rz(-2.4158106) q[2];
sx q[2];
rz(2.5570585) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.7717786) q[1];
sx q[1];
rz(-1.6057456) q[1];
sx q[1];
rz(0.82280226) q[1];
rz(-pi) q[2];
rz(0.6471031) q[3];
sx q[3];
rz(-2.3386293) q[3];
sx q[3];
rz(-0.76319198) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.42489147) q[2];
sx q[2];
rz(-0.5867914) q[2];
sx q[2];
rz(-1.1699886) q[2];
rz(0.26015002) q[3];
sx q[3];
rz(-1.5258421) q[3];
sx q[3];
rz(2.9366176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3527356) q[0];
sx q[0];
rz(-1.4997046) q[0];
sx q[0];
rz(0.66057551) q[0];
rz(-3.0109516) q[1];
sx q[1];
rz(-2.525593) q[1];
sx q[1];
rz(0.43636838) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0562387) q[0];
sx q[0];
rz(-1.071601) q[0];
sx q[0];
rz(-0.27277314) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5533668) q[2];
sx q[2];
rz(-2.4674621) q[2];
sx q[2];
rz(0.16954409) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.4727488) q[1];
sx q[1];
rz(-1.9276697) q[1];
sx q[1];
rz(-1.0297187) q[1];
rz(-1.4722669) q[3];
sx q[3];
rz(-1.2312345) q[3];
sx q[3];
rz(1.3738541) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.067430647) q[2];
sx q[2];
rz(-0.51965153) q[2];
sx q[2];
rz(2.0036428) q[2];
rz(2.2940995) q[3];
sx q[3];
rz(-1.3769826) q[3];
sx q[3];
rz(-1.725089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6013019) q[0];
sx q[0];
rz(-1.2492981) q[0];
sx q[0];
rz(-0.20703319) q[0];
rz(1.8197458) q[1];
sx q[1];
rz(-1.1489979) q[1];
sx q[1];
rz(2.102899) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10932163) q[0];
sx q[0];
rz(-2.3242852) q[0];
sx q[0];
rz(-0.82839806) q[0];
x q[1];
rz(2.8530099) q[2];
sx q[2];
rz(-0.64060694) q[2];
sx q[2];
rz(0.45329061) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.6398479) q[1];
sx q[1];
rz(-1.5540489) q[1];
sx q[1];
rz(2.5052951) q[1];
rz(-0.33090277) q[3];
sx q[3];
rz(-0.32308137) q[3];
sx q[3];
rz(-2.4289055) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.1553104) q[2];
sx q[2];
rz(-2.0955413) q[2];
sx q[2];
rz(-2.5347064) q[2];
rz(-2.8066929) q[3];
sx q[3];
rz(-1.5523942) q[3];
sx q[3];
rz(1.5437532) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63810054) q[0];
sx q[0];
rz(-1.426921) q[0];
sx q[0];
rz(-2.959429) q[0];
rz(-0.69193524) q[1];
sx q[1];
rz(-1.4357932) q[1];
sx q[1];
rz(-1.8686132) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70304027) q[0];
sx q[0];
rz(-1.5219757) q[0];
sx q[0];
rz(0.010558616) q[0];
x q[1];
rz(-2.9115306) q[2];
sx q[2];
rz(-1.4289235) q[2];
sx q[2];
rz(-1.0070247) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.0611736) q[1];
sx q[1];
rz(-0.81810942) q[1];
sx q[1];
rz(-1.2123176) q[1];
x q[2];
rz(-0.35785488) q[3];
sx q[3];
rz(-1.5113514) q[3];
sx q[3];
rz(-2.418973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.89340297) q[2];
sx q[2];
rz(-2.4290163) q[2];
sx q[2];
rz(-2.2656061) q[2];
rz(-2.857483) q[3];
sx q[3];
rz(-2.1967389) q[3];
sx q[3];
rz(-0.42377728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9559602) q[0];
sx q[0];
rz(-0.9469339) q[0];
sx q[0];
rz(-2.60485) q[0];
rz(2.5770889) q[1];
sx q[1];
rz(-0.91465488) q[1];
sx q[1];
rz(-1.5980665) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15030433) q[0];
sx q[0];
rz(-1.8841916) q[0];
sx q[0];
rz(-1.4652976) q[0];
rz(-pi) q[1];
rz(-1.364643) q[2];
sx q[2];
rz(-1.1627623) q[2];
sx q[2];
rz(-1.598151) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.1607223) q[1];
sx q[1];
rz(-1.1127143) q[1];
sx q[1];
rz(2.0715817) q[1];
rz(0.95682896) q[3];
sx q[3];
rz(-1.162863) q[3];
sx q[3];
rz(-1.8338501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.2619065) q[2];
sx q[2];
rz(-1.3085111) q[2];
sx q[2];
rz(1.7730664) q[2];
rz(2.8241099) q[3];
sx q[3];
rz(-1.2513132) q[3];
sx q[3];
rz(-0.69664636) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19318652) q[0];
sx q[0];
rz(-0.53006154) q[0];
sx q[0];
rz(-2.1004706) q[0];
rz(0.55024534) q[1];
sx q[1];
rz(-0.3717652) q[1];
sx q[1];
rz(-2.9464029) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7025906) q[0];
sx q[0];
rz(-2.9948108) q[0];
sx q[0];
rz(-0.37555666) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3615029) q[2];
sx q[2];
rz(-0.49041623) q[2];
sx q[2];
rz(1.121466) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.4976552) q[1];
sx q[1];
rz(-1.2841088) q[1];
sx q[1];
rz(3.1189256) q[1];
rz(-pi) q[2];
rz(2.6879999) q[3];
sx q[3];
rz(-2.2633584) q[3];
sx q[3];
rz(0.778824) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.13888415) q[2];
sx q[2];
rz(-1.1382853) q[2];
sx q[2];
rz(-1.2646593) q[2];
rz(1.5539315) q[3];
sx q[3];
rz(-1.9374266) q[3];
sx q[3];
rz(1.1519661) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78466648) q[0];
sx q[0];
rz(-0.21360989) q[0];
sx q[0];
rz(-2.7333562) q[0];
rz(-1.6881855) q[1];
sx q[1];
rz(-1.6994349) q[1];
sx q[1];
rz(2.2834987) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8342339) q[0];
sx q[0];
rz(-1.8862794) q[0];
sx q[0];
rz(0.90641109) q[0];
rz(-pi) q[1];
rz(-1.6995605) q[2];
sx q[2];
rz(-1.4742645) q[2];
sx q[2];
rz(1.0573204) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.1979023) q[1];
sx q[1];
rz(-2.5630782) q[1];
sx q[1];
rz(0.79165081) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.016641) q[3];
sx q[3];
rz(-0.88192421) q[3];
sx q[3];
rz(0.75923789) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.475829) q[2];
sx q[2];
rz(-0.88226157) q[2];
sx q[2];
rz(2.8909454) q[2];
rz(2.2541239) q[3];
sx q[3];
rz(-1.1424501) q[3];
sx q[3];
rz(0.34621507) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86824054) q[0];
sx q[0];
rz(-0.98889416) q[0];
sx q[0];
rz(0.86531472) q[0];
rz(-2.6189651) q[1];
sx q[1];
rz(-2.3326383) q[1];
sx q[1];
rz(1.685198) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9495372) q[0];
sx q[0];
rz(-2.9350781) q[0];
sx q[0];
rz(-1.2386991) q[0];
rz(-pi) q[1];
rz(-0.9002879) q[2];
sx q[2];
rz(-1.4354726) q[2];
sx q[2];
rz(1.6953261) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.8700614) q[1];
sx q[1];
rz(-1.2200331) q[1];
sx q[1];
rz(-0.49141541) q[1];
x q[2];
rz(0.51442149) q[3];
sx q[3];
rz(-1.2687917) q[3];
sx q[3];
rz(2.746904) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.3978079) q[2];
sx q[2];
rz(-0.91181552) q[2];
sx q[2];
rz(1.816681) q[2];
rz(-1.2815453) q[3];
sx q[3];
rz(-2.2478588) q[3];
sx q[3];
rz(-0.76461422) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21541967) q[0];
sx q[0];
rz(-1.3289435) q[0];
sx q[0];
rz(-1.5935612) q[0];
rz(-0.46319766) q[1];
sx q[1];
rz(-1.5402272) q[1];
sx q[1];
rz(2.872749) q[1];
rz(1.8824957) q[2];
sx q[2];
rz(-2.1915956) q[2];
sx q[2];
rz(-1.6804463) q[2];
rz(-2.4687121) q[3];
sx q[3];
rz(-1.665848) q[3];
sx q[3];
rz(3.0116826) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
