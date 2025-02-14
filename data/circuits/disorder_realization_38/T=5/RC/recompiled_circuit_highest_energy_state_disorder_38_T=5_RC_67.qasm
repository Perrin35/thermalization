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
rz(0.44904798) q[0];
sx q[0];
rz(-2.0023876) q[0];
sx q[0];
rz(0.44266242) q[0];
rz(-2.9709385) q[1];
sx q[1];
rz(-0.79167384) q[1];
sx q[1];
rz(-0.72905529) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3826844) q[0];
sx q[0];
rz(-2.2331862) q[0];
sx q[0];
rz(1.054396) q[0];
x q[1];
rz(-1.5698956) q[2];
sx q[2];
rz(-1.1810762) q[2];
sx q[2];
rz(0.87504059) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.52965611) q[1];
sx q[1];
rz(-1.3289598) q[1];
sx q[1];
rz(0.87301689) q[1];
rz(-1.0341845) q[3];
sx q[3];
rz(-0.70176673) q[3];
sx q[3];
rz(-0.90493363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.2810104) q[2];
sx q[2];
rz(-0.082110114) q[2];
sx q[2];
rz(2.3561467) q[2];
rz(2.1752518) q[3];
sx q[3];
rz(-2.0735425) q[3];
sx q[3];
rz(2.5016224) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6099089) q[0];
sx q[0];
rz(-2.5737679) q[0];
sx q[0];
rz(0.63572836) q[0];
rz(0.6334148) q[1];
sx q[1];
rz(-2.3737291) q[1];
sx q[1];
rz(-2.8107218) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2658871) q[0];
sx q[0];
rz(-0.71015753) q[0];
sx q[0];
rz(0.52010923) q[0];
rz(-1.517968) q[2];
sx q[2];
rz(-0.46762782) q[2];
sx q[2];
rz(0.29948452) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0111661) q[1];
sx q[1];
rz(-3.0741709) q[1];
sx q[1];
rz(-1.0481337) q[1];
rz(-pi) q[2];
x q[2];
rz(0.084609065) q[3];
sx q[3];
rz(-1.5584018) q[3];
sx q[3];
rz(0.4735213) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.6364381) q[2];
sx q[2];
rz(-0.81710368) q[2];
sx q[2];
rz(-0.46112296) q[2];
rz(2.4165706) q[3];
sx q[3];
rz(-0.45261639) q[3];
sx q[3];
rz(0.66864526) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.649699) q[0];
sx q[0];
rz(-1.3839586) q[0];
sx q[0];
rz(-2.6032676) q[0];
rz(0.0093983924) q[1];
sx q[1];
rz(-2.5098269) q[1];
sx q[1];
rz(-1.9422772) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6865348) q[0];
sx q[0];
rz(-1.063534) q[0];
sx q[0];
rz(0.051183609) q[0];
rz(-pi) q[1];
rz(-0.17260562) q[2];
sx q[2];
rz(-2.4510018) q[2];
sx q[2];
rz(2.0334854) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.6965683) q[1];
sx q[1];
rz(-0.68587947) q[1];
sx q[1];
rz(2.1522151) q[1];
rz(-2.2975397) q[3];
sx q[3];
rz(-1.7622602) q[3];
sx q[3];
rz(1.307631) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.37375307) q[2];
sx q[2];
rz(-1.6497296) q[2];
sx q[2];
rz(-0.72244942) q[2];
rz(-2.5433698) q[3];
sx q[3];
rz(-0.0349508) q[3];
sx q[3];
rz(1.9237062) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5825321) q[0];
sx q[0];
rz(-1.2578955) q[0];
sx q[0];
rz(-0.019056888) q[0];
rz(-1.6825153) q[1];
sx q[1];
rz(-0.2225114) q[1];
sx q[1];
rz(-0.51024514) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.048240926) q[0];
sx q[0];
rz(-2.4722795) q[0];
sx q[0];
rz(1.1360243) q[0];
rz(-pi) q[1];
rz(-2.7925582) q[2];
sx q[2];
rz(-2.3885997) q[2];
sx q[2];
rz(-0.53607625) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.45106798) q[1];
sx q[1];
rz(-2.1803586) q[1];
sx q[1];
rz(-1.068348) q[1];
rz(-pi) q[2];
rz(0.31655772) q[3];
sx q[3];
rz(-0.43310851) q[3];
sx q[3];
rz(0.18726997) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9973008) q[2];
sx q[2];
rz(-1.68953) q[2];
sx q[2];
rz(-2.9193381) q[2];
rz(0.079515919) q[3];
sx q[3];
rz(-2.5955718) q[3];
sx q[3];
rz(-2.5047746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52614373) q[0];
sx q[0];
rz(-1.0960824) q[0];
sx q[0];
rz(1.3413606) q[0];
rz(0.51097956) q[1];
sx q[1];
rz(-1.363441) q[1];
sx q[1];
rz(2.708639) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.88192542) q[0];
sx q[0];
rz(-2.8779463) q[0];
sx q[0];
rz(2.9211098) q[0];
rz(-2.9278361) q[2];
sx q[2];
rz(-1.1756437) q[2];
sx q[2];
rz(-1.9825801) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.67636) q[1];
sx q[1];
rz(-1.3238878) q[1];
sx q[1];
rz(-1.9104772) q[1];
rz(2.8540887) q[3];
sx q[3];
rz(-2.3477481) q[3];
sx q[3];
rz(2.7903729) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0114835) q[2];
sx q[2];
rz(-2.2138962) q[2];
sx q[2];
rz(-2.1437272) q[2];
rz(-2.4326371) q[3];
sx q[3];
rz(-0.38781375) q[3];
sx q[3];
rz(-0.97878218) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3157432) q[0];
sx q[0];
rz(-2.5542927) q[0];
sx q[0];
rz(-0.89383268) q[0];
rz(2.1122487) q[1];
sx q[1];
rz(-0.7411595) q[1];
sx q[1];
rz(-3.0267874) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.288088) q[0];
sx q[0];
rz(-1.813632) q[0];
sx q[0];
rz(-0.22260861) q[0];
x q[1];
rz(0.64115142) q[2];
sx q[2];
rz(-1.8617407) q[2];
sx q[2];
rz(3.0464622) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.5789168) q[1];
sx q[1];
rz(-1.8496184) q[1];
sx q[1];
rz(-1.2827966) q[1];
rz(-1.6345665) q[3];
sx q[3];
rz(-1.8511417) q[3];
sx q[3];
rz(0.75322039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.7041695) q[2];
sx q[2];
rz(-0.88714522) q[2];
sx q[2];
rz(-0.21370055) q[2];
rz(-2.5050488) q[3];
sx q[3];
rz(-2.611219) q[3];
sx q[3];
rz(0.73591939) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32695025) q[0];
sx q[0];
rz(-0.7433759) q[0];
sx q[0];
rz(-3.0315234) q[0];
rz(-3.0649109) q[1];
sx q[1];
rz(-2.2814467) q[1];
sx q[1];
rz(2.0032047) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6784043) q[0];
sx q[0];
rz(-0.41751719) q[0];
sx q[0];
rz(-0.56708401) q[0];
rz(-pi) q[1];
rz(-0.0037527379) q[2];
sx q[2];
rz(-0.787799) q[2];
sx q[2];
rz(-1.4472199) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.0323022) q[1];
sx q[1];
rz(-1.7545953) q[1];
sx q[1];
rz(-1.079783) q[1];
x q[2];
rz(1.3567233) q[3];
sx q[3];
rz(-1.8118637) q[3];
sx q[3];
rz(1.1837219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.25973213) q[2];
sx q[2];
rz(-2.0596108) q[2];
sx q[2];
rz(-2.760375) q[2];
rz(3.0242053) q[3];
sx q[3];
rz(-2.6365247) q[3];
sx q[3];
rz(0.086732619) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36670244) q[0];
sx q[0];
rz(-2.369304) q[0];
sx q[0];
rz(0.49852398) q[0];
rz(-2.3122834) q[1];
sx q[1];
rz(-2.3259951) q[1];
sx q[1];
rz(0.097749762) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5437336) q[0];
sx q[0];
rz(-1.8761138) q[0];
sx q[0];
rz(-2.0241364) q[0];
rz(-2.8016053) q[2];
sx q[2];
rz(-1.7113216) q[2];
sx q[2];
rz(-0.7896151) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.5796639) q[1];
sx q[1];
rz(-1.0439928) q[1];
sx q[1];
rz(-2.4234523) q[1];
rz(-pi) q[2];
rz(1.5427715) q[3];
sx q[3];
rz(-0.32169473) q[3];
sx q[3];
rz(3.0672586) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.8349614) q[2];
sx q[2];
rz(-1.9913048) q[2];
sx q[2];
rz(0.59070307) q[2];
rz(-3.0513884) q[3];
sx q[3];
rz(-2.7019751) q[3];
sx q[3];
rz(2.2265767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0331994) q[0];
sx q[0];
rz(-2.2985701) q[0];
sx q[0];
rz(2.3868308) q[0];
rz(2.8056878) q[1];
sx q[1];
rz(-0.55481189) q[1];
sx q[1];
rz(1.0364484) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91387902) q[0];
sx q[0];
rz(-0.69509414) q[0];
sx q[0];
rz(0.18385017) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.039128379) q[2];
sx q[2];
rz(-0.48155287) q[2];
sx q[2];
rz(2.0619543) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.4200585) q[1];
sx q[1];
rz(-2.0871401) q[1];
sx q[1];
rz(1.615585) q[1];
rz(1.8004981) q[3];
sx q[3];
rz(-1.5840845) q[3];
sx q[3];
rz(-2.5984327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-3.0466517) q[2];
sx q[2];
rz(-1.4163821) q[2];
sx q[2];
rz(2.6851192) q[2];
rz(0.60232317) q[3];
sx q[3];
rz(-2.9538302) q[3];
sx q[3];
rz(-2.203234) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25191054) q[0];
sx q[0];
rz(-1.849181) q[0];
sx q[0];
rz(0.46345261) q[0];
rz(3.1262596) q[1];
sx q[1];
rz(-0.066611193) q[1];
sx q[1];
rz(-2.8212246) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25523089) q[0];
sx q[0];
rz(-1.4275274) q[0];
sx q[0];
rz(-1.6485639) q[0];
x q[1];
rz(-2.4762829) q[2];
sx q[2];
rz(-2.7501366) q[2];
sx q[2];
rz(-0.91400439) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5825504) q[1];
sx q[1];
rz(-2.1782785) q[1];
sx q[1];
rz(2.5431125) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4418751) q[3];
sx q[3];
rz(-2.7072002) q[3];
sx q[3];
rz(0.35067973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(3.1128803) q[2];
sx q[2];
rz(-2.4148603) q[2];
sx q[2];
rz(-1.0090562) q[2];
rz(-2.3446337) q[3];
sx q[3];
rz(-1.2538486) q[3];
sx q[3];
rz(2.8033281) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0924031) q[0];
sx q[0];
rz(-1.6272463) q[0];
sx q[0];
rz(1.8820681) q[0];
rz(0.084820329) q[1];
sx q[1];
rz(-1.6592574) q[1];
sx q[1];
rz(1.4316373) q[1];
rz(1.3289497) q[2];
sx q[2];
rz(-1.5089265) q[2];
sx q[2];
rz(-0.93969719) q[2];
rz(2.6081035) q[3];
sx q[3];
rz(-0.34779741) q[3];
sx q[3];
rz(1.6602914) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
