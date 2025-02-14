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
rz(0.64492172) q[0];
sx q[0];
rz(2.6725197) q[0];
sx q[0];
rz(7.2735431) q[0];
rz(-1.1733836) q[1];
sx q[1];
rz(2.2819509) q[1];
sx q[1];
rz(11.717164) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3135017) q[0];
sx q[0];
rz(-1.573631) q[0];
sx q[0];
rz(1.6046197) q[0];
x q[1];
rz(1.1728549) q[2];
sx q[2];
rz(-2.401377) q[2];
sx q[2];
rz(0.55392439) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7265004) q[1];
sx q[1];
rz(-2.1884568) q[1];
sx q[1];
rz(-2.8128977) q[1];
x q[2];
rz(-1.7683214) q[3];
sx q[3];
rz(-0.53086262) q[3];
sx q[3];
rz(-1.498675) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.2363756) q[2];
sx q[2];
rz(-1.6340916) q[2];
sx q[2];
rz(-0.53984731) q[2];
rz(-1.5652462) q[3];
sx q[3];
rz(-0.40894517) q[3];
sx q[3];
rz(1.7319771) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
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
rz(-0.09963116) q[0];
sx q[0];
rz(-2.4653682) q[0];
sx q[0];
rz(-1.8145632) q[0];
rz(-0.78500336) q[1];
sx q[1];
rz(-1.7186586) q[1];
sx q[1];
rz(2.6410417) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1091565) q[0];
sx q[0];
rz(-1.5958428) q[0];
sx q[0];
rz(-0.5513726) q[0];
rz(-1.9376061) q[2];
sx q[2];
rz(-1.3932835) q[2];
sx q[2];
rz(0.29581636) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(3.0620898) q[1];
sx q[1];
rz(-0.70598733) q[1];
sx q[1];
rz(1.0393618) q[1];
rz(-0.80970069) q[3];
sx q[3];
rz(-1.2375583) q[3];
sx q[3];
rz(2.7908195) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.0850204) q[2];
sx q[2];
rz(-2.4041921) q[2];
sx q[2];
rz(2.6596587) q[2];
rz(0.36007544) q[3];
sx q[3];
rz(-1.1432546) q[3];
sx q[3];
rz(2.9329407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8423186) q[0];
sx q[0];
rz(-1.1085008) q[0];
sx q[0];
rz(-0.47766787) q[0];
rz(-0.85992366) q[1];
sx q[1];
rz(-2.0932902) q[1];
sx q[1];
rz(-0.28712505) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56518131) q[0];
sx q[0];
rz(-2.9991096) q[0];
sx q[0];
rz(-1.0866685) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7870891) q[2];
sx q[2];
rz(-1.7743127) q[2];
sx q[2];
rz(3.1151724) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.7309235) q[1];
sx q[1];
rz(-0.6836709) q[1];
sx q[1];
rz(1.5846496) q[1];
rz(-pi) q[2];
x q[2];
rz(0.085137376) q[3];
sx q[3];
rz(-1.9508024) q[3];
sx q[3];
rz(0.054084965) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.75256789) q[2];
sx q[2];
rz(-1.7914881) q[2];
sx q[2];
rz(3.1217421) q[2];
rz(-2.1974473) q[3];
sx q[3];
rz(-2.4343334) q[3];
sx q[3];
rz(2.7437362) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95553628) q[0];
sx q[0];
rz(-2.5113386) q[0];
sx q[0];
rz(-1.3007042) q[0];
rz(-2.6553254) q[1];
sx q[1];
rz(-1.174289) q[1];
sx q[1];
rz(0.035331443) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.080950532) q[0];
sx q[0];
rz(-1.5968906) q[0];
sx q[0];
rz(1.1423201) q[0];
rz(-pi) q[1];
rz(-2.6258518) q[2];
sx q[2];
rz(-0.44937849) q[2];
sx q[2];
rz(0.55381227) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.3934734) q[1];
sx q[1];
rz(-1.5474209) q[1];
sx q[1];
rz(1.5651902) q[1];
x q[2];
rz(0.62832812) q[3];
sx q[3];
rz(-1.105996) q[3];
sx q[3];
rz(-0.94098722) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.6639634) q[2];
sx q[2];
rz(-2.5783381) q[2];
sx q[2];
rz(-0.25838724) q[2];
rz(2.431331) q[3];
sx q[3];
rz(-0.60451549) q[3];
sx q[3];
rz(1.5593504) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48359394) q[0];
sx q[0];
rz(-0.76376629) q[0];
sx q[0];
rz(-0.64436954) q[0];
rz(-2.1517892) q[1];
sx q[1];
rz(-0.80392307) q[1];
sx q[1];
rz(-2.03233) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4424446) q[0];
sx q[0];
rz(-2.0936396) q[0];
sx q[0];
rz(0.22626489) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.039626683) q[2];
sx q[2];
rz(-1.0575231) q[2];
sx q[2];
rz(0.64943343) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.4586045) q[1];
sx q[1];
rz(-2.3420077) q[1];
sx q[1];
rz(1.7728642) q[1];
rz(-pi) q[2];
rz(2.3886834) q[3];
sx q[3];
rz(-1.3594106) q[3];
sx q[3];
rz(1.961238) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.41205078) q[2];
sx q[2];
rz(-1.5770301) q[2];
sx q[2];
rz(0.61005074) q[2];
rz(-0.080502056) q[3];
sx q[3];
rz(-0.1463612) q[3];
sx q[3];
rz(-0.0065053594) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7215111) q[0];
sx q[0];
rz(-0.93669909) q[0];
sx q[0];
rz(-1.4790081) q[0];
rz(1.9526019) q[1];
sx q[1];
rz(-2.7956796) q[1];
sx q[1];
rz(-1.6993274) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8540031) q[0];
sx q[0];
rz(-2.1882399) q[0];
sx q[0];
rz(-2.94728) q[0];
rz(-0.58143388) q[2];
sx q[2];
rz(-1.6260885) q[2];
sx q[2];
rz(-0.96382574) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.4870781) q[1];
sx q[1];
rz(-2.2527825) q[1];
sx q[1];
rz(-0.42907663) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9154432) q[3];
sx q[3];
rz(-2.0807869) q[3];
sx q[3];
rz(0.6450212) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.46001616) q[2];
sx q[2];
rz(-0.67395335) q[2];
sx q[2];
rz(-2.4665534) q[2];
rz(0.33440822) q[3];
sx q[3];
rz(-1.6655917) q[3];
sx q[3];
rz(-2.7736751) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6595031) q[0];
sx q[0];
rz(-0.98564321) q[0];
sx q[0];
rz(-0.12568812) q[0];
rz(0.19371678) q[1];
sx q[1];
rz(-0.88231641) q[1];
sx q[1];
rz(1.151459) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4070258) q[0];
sx q[0];
rz(-2.1687963) q[0];
sx q[0];
rz(2.1728188) q[0];
rz(-pi) q[1];
rz(2.861768) q[2];
sx q[2];
rz(-1.9283659) q[2];
sx q[2];
rz(1.2841061) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9509778) q[1];
sx q[1];
rz(-2.0101317) q[1];
sx q[1];
rz(2.5530035) q[1];
x q[2];
rz(1.098387) q[3];
sx q[3];
rz(-2.189872) q[3];
sx q[3];
rz(0.38451871) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.2842399) q[2];
sx q[2];
rz(-1.7809296) q[2];
sx q[2];
rz(-0.24687684) q[2];
rz(-0.85865584) q[3];
sx q[3];
rz(-0.98444986) q[3];
sx q[3];
rz(-2.8624559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3566647) q[0];
sx q[0];
rz(-0.44342884) q[0];
sx q[0];
rz(-0.32522935) q[0];
rz(-0.52109703) q[1];
sx q[1];
rz(-0.63322133) q[1];
sx q[1];
rz(2.8281143) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3081261) q[0];
sx q[0];
rz(-2.0692056) q[0];
sx q[0];
rz(-0.68584401) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3455079) q[2];
sx q[2];
rz(-1.5048001) q[2];
sx q[2];
rz(-2.5661732) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.5446452) q[1];
sx q[1];
rz(-1.2352422) q[1];
sx q[1];
rz(-0.38864522) q[1];
rz(-pi) q[2];
x q[2];
rz(0.20200396) q[3];
sx q[3];
rz(-1.0583377) q[3];
sx q[3];
rz(-0.92512586) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.7079805) q[2];
sx q[2];
rz(-1.2182451) q[2];
sx q[2];
rz(1.8617967) q[2];
rz(2.4158939) q[3];
sx q[3];
rz(-1.1596707) q[3];
sx q[3];
rz(2.3316135) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5028266) q[0];
sx q[0];
rz(-1.9483197) q[0];
sx q[0];
rz(0.1499114) q[0];
rz(1.2431078) q[1];
sx q[1];
rz(-1.5877692) q[1];
sx q[1];
rz(-1.6601723) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.178326) q[0];
sx q[0];
rz(-1.5218922) q[0];
sx q[0];
rz(0.045128926) q[0];
x q[1];
rz(-1.7779839) q[2];
sx q[2];
rz(-2.8958791) q[2];
sx q[2];
rz(2.4893275) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.31563035) q[1];
sx q[1];
rz(-0.44188979) q[1];
sx q[1];
rz(0.52200861) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0143901) q[3];
sx q[3];
rz(-1.8133928) q[3];
sx q[3];
rz(-0.85240817) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.82568613) q[2];
sx q[2];
rz(-1.3797727) q[2];
sx q[2];
rz(2.1659577) q[2];
rz(2.0129096) q[3];
sx q[3];
rz(-0.55207878) q[3];
sx q[3];
rz(-0.25477195) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16633701) q[0];
sx q[0];
rz(-2.14125) q[0];
sx q[0];
rz(-2.9794203) q[0];
rz(1.3316679) q[1];
sx q[1];
rz(-0.63260308) q[1];
sx q[1];
rz(0.046028927) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9105658) q[0];
sx q[0];
rz(-1.5648769) q[0];
sx q[0];
rz(1.5886515) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9580905) q[2];
sx q[2];
rz(-2.0047054) q[2];
sx q[2];
rz(0.48112416) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-3.0303453) q[1];
sx q[1];
rz(-1.971389) q[1];
sx q[1];
rz(2.9811893) q[1];
rz(-pi) q[2];
rz(1.7415984) q[3];
sx q[3];
rz(-1.2360337) q[3];
sx q[3];
rz(0.46250175) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.33596805) q[2];
sx q[2];
rz(-0.38463548) q[2];
sx q[2];
rz(1.1495205) q[2];
rz(0.51291054) q[3];
sx q[3];
rz(-1.7549113) q[3];
sx q[3];
rz(-2.8368867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.649986) q[0];
sx q[0];
rz(-2.2254324) q[0];
sx q[0];
rz(-1.9944763) q[0];
rz(-0.06123771) q[1];
sx q[1];
rz(-1.1920659) q[1];
sx q[1];
rz(1.7658284) q[1];
rz(-2.1381151) q[2];
sx q[2];
rz(-0.85366953) q[2];
sx q[2];
rz(0.015346957) q[2];
rz(0.41027222) q[3];
sx q[3];
rz(-0.82220746) q[3];
sx q[3];
rz(-2.0668277) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
