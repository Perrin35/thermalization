OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.2849046) q[0];
sx q[0];
rz(-1.8104799) q[0];
sx q[0];
rz(2.3321505) q[0];
rz(0.59960214) q[1];
sx q[1];
rz(-2.0166346) q[1];
sx q[1];
rz(-0.52176276) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7358462) q[0];
sx q[0];
rz(-0.1400811) q[0];
sx q[0];
rz(1.5790042) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2748364) q[2];
sx q[2];
rz(-2.1115536) q[2];
sx q[2];
rz(2.7014033) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.45568902) q[1];
sx q[1];
rz(-2.0200811) q[1];
sx q[1];
rz(2.8390769) q[1];
rz(-pi) q[2];
rz(-1.4843602) q[3];
sx q[3];
rz(-2.1051791) q[3];
sx q[3];
rz(-0.57356452) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.9894422) q[2];
sx q[2];
rz(-0.70789727) q[2];
sx q[2];
rz(1.7789827) q[2];
rz(1.4710434) q[3];
sx q[3];
rz(-1.5444642) q[3];
sx q[3];
rz(0.41489261) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3926369) q[0];
sx q[0];
rz(-1.7935268) q[0];
sx q[0];
rz(1.0697399) q[0];
rz(2.3906129) q[1];
sx q[1];
rz(-0.90194482) q[1];
sx q[1];
rz(0.78838563) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8892944) q[0];
sx q[0];
rz(-1.3836593) q[0];
sx q[0];
rz(-0.62690134) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4761202) q[2];
sx q[2];
rz(-1.62684) q[2];
sx q[2];
rz(-1.9201345) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.3020116) q[1];
sx q[1];
rz(-0.55802781) q[1];
sx q[1];
rz(1.9424979) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0072924) q[3];
sx q[3];
rz(-1.1009163) q[3];
sx q[3];
rz(1.7738455) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.02701935) q[2];
sx q[2];
rz(-2.7832649) q[2];
sx q[2];
rz(2.2349854) q[2];
rz(1.00057) q[3];
sx q[3];
rz(-1.118719) q[3];
sx q[3];
rz(-0.51893836) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24666102) q[0];
sx q[0];
rz(-2.7592359) q[0];
sx q[0];
rz(-1.83778) q[0];
rz(1.9815824) q[1];
sx q[1];
rz(-1.8142533) q[1];
sx q[1];
rz(-1.4132285) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.055186836) q[0];
sx q[0];
rz(-0.16567437) q[0];
sx q[0];
rz(-0.55060668) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.5858461) q[2];
sx q[2];
rz(-0.72524643) q[2];
sx q[2];
rz(0.21373978) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.47067898) q[1];
sx q[1];
rz(-1.7470198) q[1];
sx q[1];
rz(-2.6754968) q[1];
rz(-pi) q[2];
x q[2];
rz(1.71536) q[3];
sx q[3];
rz(-0.72025296) q[3];
sx q[3];
rz(-1.6584058) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.7580737) q[2];
sx q[2];
rz(-1.588409) q[2];
sx q[2];
rz(0.12913945) q[2];
rz(0.50390759) q[3];
sx q[3];
rz(-1.7241071) q[3];
sx q[3];
rz(-2.5655139) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.031216089) q[0];
sx q[0];
rz(-1.2126558) q[0];
sx q[0];
rz(-0.84621286) q[0];
rz(-0.57202488) q[1];
sx q[1];
rz(-1.5049479) q[1];
sx q[1];
rz(-1.2342854) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18674984) q[0];
sx q[0];
rz(-0.34663793) q[0];
sx q[0];
rz(2.8524272) q[0];
rz(-pi) q[1];
rz(-2.8568966) q[2];
sx q[2];
rz(-2.0139593) q[2];
sx q[2];
rz(-2.5131132) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.438414) q[1];
sx q[1];
rz(-1.6513737) q[1];
sx q[1];
rz(2.9318649) q[1];
rz(-2.0999072) q[3];
sx q[3];
rz(-2.2941781) q[3];
sx q[3];
rz(0.46284562) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7531551) q[2];
sx q[2];
rz(-2.2074102) q[2];
sx q[2];
rz(1.4446806) q[2];
rz(1.2706903) q[3];
sx q[3];
rz(-1.5710187) q[3];
sx q[3];
rz(1.1893893) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5020849) q[0];
sx q[0];
rz(-1.4690228) q[0];
sx q[0];
rz(-2.0462659) q[0];
rz(3.0785839) q[1];
sx q[1];
rz(-1.4728225) q[1];
sx q[1];
rz(-2.1953886) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99710714) q[0];
sx q[0];
rz(-1.5427663) q[0];
sx q[0];
rz(1.0260242) q[0];
rz(-pi) q[1];
rz(1.9502207) q[2];
sx q[2];
rz(-2.1721031) q[2];
sx q[2];
rz(1.7446818) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.9506388) q[1];
sx q[1];
rz(-2.6304448) q[1];
sx q[1];
rz(0.010015476) q[1];
rz(-pi) q[2];
rz(1.9619476) q[3];
sx q[3];
rz(-1.3538651) q[3];
sx q[3];
rz(-0.85238487) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.63518628) q[2];
sx q[2];
rz(-1.0794285) q[2];
sx q[2];
rz(1.5080473) q[2];
rz(-0.15277282) q[3];
sx q[3];
rz(-1.907405) q[3];
sx q[3];
rz(-0.93301409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2372811) q[0];
sx q[0];
rz(-0.58552423) q[0];
sx q[0];
rz(-0.10865077) q[0];
rz(-3.0142504) q[1];
sx q[1];
rz(-1.8904949) q[1];
sx q[1];
rz(2.2551575) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21289794) q[0];
sx q[0];
rz(-1.587364) q[0];
sx q[0];
rz(2.6221656) q[0];
rz(-pi) q[1];
rz(-1.3662874) q[2];
sx q[2];
rz(-1.5123774) q[2];
sx q[2];
rz(-1.2893334) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.8005714) q[1];
sx q[1];
rz(-1.6280975) q[1];
sx q[1];
rz(-1.351746) q[1];
rz(0.57678595) q[3];
sx q[3];
rz(-2.9401719) q[3];
sx q[3];
rz(-0.52514986) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.2966557) q[2];
sx q[2];
rz(-0.79534328) q[2];
sx q[2];
rz(0.33357683) q[2];
rz(-0.9345471) q[3];
sx q[3];
rz(-2.3881674) q[3];
sx q[3];
rz(2.2540588) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4535256) q[0];
sx q[0];
rz(-1.6418566) q[0];
sx q[0];
rz(-2.8330579) q[0];
rz(2.535179) q[1];
sx q[1];
rz(-2.2589222) q[1];
sx q[1];
rz(0.44953129) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.318293) q[0];
sx q[0];
rz(-1.6345981) q[0];
sx q[0];
rz(-2.6266111) q[0];
rz(2.6753898) q[2];
sx q[2];
rz(-0.84377366) q[2];
sx q[2];
rz(-0.94949338) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.3492185) q[1];
sx q[1];
rz(-1.4687531) q[1];
sx q[1];
rz(2.1898063) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.94664871) q[3];
sx q[3];
rz(-2.3376541) q[3];
sx q[3];
rz(-2.338666) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.061607925) q[2];
sx q[2];
rz(-1.6500762) q[2];
sx q[2];
rz(-1.1770581) q[2];
rz(2.4401149) q[3];
sx q[3];
rz(-0.25291118) q[3];
sx q[3];
rz(-2.285932) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6845282) q[0];
sx q[0];
rz(-0.51790154) q[0];
sx q[0];
rz(1.49217) q[0];
rz(-1.0555142) q[1];
sx q[1];
rz(-1.5558473) q[1];
sx q[1];
rz(-0.42748705) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1221177) q[0];
sx q[0];
rz(-0.97163288) q[0];
sx q[0];
rz(0.63728665) q[0];
rz(-pi) q[1];
rz(-0.43085499) q[2];
sx q[2];
rz(-1.4240032) q[2];
sx q[2];
rz(-2.1374747) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.44459773) q[1];
sx q[1];
rz(-2.8279404) q[1];
sx q[1];
rz(0.072254953) q[1];
x q[2];
rz(-1.2990789) q[3];
sx q[3];
rz(-1.875056) q[3];
sx q[3];
rz(0.4189241) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.9588354) q[2];
sx q[2];
rz(-1.8104825) q[2];
sx q[2];
rz(-1.9367564) q[2];
rz(3.0086573) q[3];
sx q[3];
rz(-3.0429621) q[3];
sx q[3];
rz(-2.0601864) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28250113) q[0];
sx q[0];
rz(-1.7318672) q[0];
sx q[0];
rz(0.45243636) q[0];
rz(-2.2826507) q[1];
sx q[1];
rz(-2.5622538) q[1];
sx q[1];
rz(-1.6890242) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5256663) q[0];
sx q[0];
rz(-2.1480664) q[0];
sx q[0];
rz(0.99212118) q[0];
x q[1];
rz(0.33272393) q[2];
sx q[2];
rz(-2.1004268) q[2];
sx q[2];
rz(0.67654787) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(3.0124546) q[1];
sx q[1];
rz(-2.635639) q[1];
sx q[1];
rz(1.5652191) q[1];
x q[2];
rz(2.3708569) q[3];
sx q[3];
rz(-0.49388921) q[3];
sx q[3];
rz(-2.1300742) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.9026044) q[2];
sx q[2];
rz(-1.8391823) q[2];
sx q[2];
rz(0.39546173) q[2];
rz(-2.0579193) q[3];
sx q[3];
rz(-2.5189221) q[3];
sx q[3];
rz(-1.4130939) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.285242) q[0];
sx q[0];
rz(-3.0265891) q[0];
sx q[0];
rz(1.2913936) q[0];
rz(2.3639288) q[1];
sx q[1];
rz(-1.5592557) q[1];
sx q[1];
rz(1.8596328) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4829452) q[0];
sx q[0];
rz(-1.7883669) q[0];
sx q[0];
rz(-1.3402142) q[0];
rz(-pi) q[1];
rz(-2.8516232) q[2];
sx q[2];
rz(-2.4735138) q[2];
sx q[2];
rz(-1.5886605) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.3471372) q[1];
sx q[1];
rz(-1.0442808) q[1];
sx q[1];
rz(-1.305425) q[1];
x q[2];
rz(-3.0324803) q[3];
sx q[3];
rz(-1.0541354) q[3];
sx q[3];
rz(-1.9369879) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.164244) q[2];
sx q[2];
rz(-2.7898596) q[2];
sx q[2];
rz(0.39883167) q[2];
rz(0.93419689) q[3];
sx q[3];
rz(-2.4103006) q[3];
sx q[3];
rz(-3.0493693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.6964523) q[0];
sx q[0];
rz(-0.91309375) q[0];
sx q[0];
rz(0.0050807411) q[0];
rz(2.483881) q[1];
sx q[1];
rz(-1.7291768) q[1];
sx q[1];
rz(-1.9531858) q[1];
rz(2.3376113) q[2];
sx q[2];
rz(-1.3406499) q[2];
sx q[2];
rz(-1.7743877) q[2];
rz(-1.4767968) q[3];
sx q[3];
rz(-0.58782676) q[3];
sx q[3];
rz(1.6997433) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
