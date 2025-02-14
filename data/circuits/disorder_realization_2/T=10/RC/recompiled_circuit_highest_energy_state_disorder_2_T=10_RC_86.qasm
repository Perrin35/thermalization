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
rz(2.958137) q[0];
sx q[0];
rz(-2.3975211) q[0];
sx q[0];
rz(-2.0896572) q[0];
rz(0.19368859) q[1];
sx q[1];
rz(2.5084578) q[1];
sx q[1];
rz(8.8517744) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.075112315) q[0];
sx q[0];
rz(-1.6329935) q[0];
sx q[0];
rz(-0.94554995) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.139976) q[2];
sx q[2];
rz(-0.9316906) q[2];
sx q[2];
rz(2.4205645) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.68202735) q[1];
sx q[1];
rz(-0.50217705) q[1];
sx q[1];
rz(-0.14634303) q[1];
rz(-pi) q[2];
rz(-1.6735733) q[3];
sx q[3];
rz(-0.87689084) q[3];
sx q[3];
rz(-0.68460195) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.8594592) q[2];
sx q[2];
rz(-2.8318475) q[2];
sx q[2];
rz(-1.0563043) q[2];
rz(-3.1193962) q[3];
sx q[3];
rz(-2.3789417) q[3];
sx q[3];
rz(1.7215884) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9098772) q[0];
sx q[0];
rz(-0.48119369) q[0];
sx q[0];
rz(2.7114482) q[0];
rz(3.0138956) q[1];
sx q[1];
rz(-1.1558497) q[1];
sx q[1];
rz(1.7040303) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6352954) q[0];
sx q[0];
rz(-2.7397635) q[0];
sx q[0];
rz(-2.4319629) q[0];
x q[1];
rz(0.050927863) q[2];
sx q[2];
rz(-1.5136592) q[2];
sx q[2];
rz(0.47540755) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.3900657) q[1];
sx q[1];
rz(-1.2428778) q[1];
sx q[1];
rz(-2.4310914) q[1];
rz(-pi) q[2];
x q[2];
rz(0.53967472) q[3];
sx q[3];
rz(-1.7575348) q[3];
sx q[3];
rz(0.58247551) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(3.0251856) q[2];
sx q[2];
rz(-1.6802639) q[2];
sx q[2];
rz(2.6961668) q[2];
rz(2.7323501) q[3];
sx q[3];
rz(-1.0990812) q[3];
sx q[3];
rz(-2.2621431) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5889848) q[0];
sx q[0];
rz(-3.0486076) q[0];
sx q[0];
rz(0.71480042) q[0];
rz(-2.0843166) q[1];
sx q[1];
rz(-0.42066586) q[1];
sx q[1];
rz(2.8143299) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3540368) q[0];
sx q[0];
rz(-1.7374514) q[0];
sx q[0];
rz(2.119407) q[0];
rz(-pi) q[1];
rz(-0.94560854) q[2];
sx q[2];
rz(-0.66788061) q[2];
sx q[2];
rz(2.5651074) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.2043874) q[1];
sx q[1];
rz(-1.4193077) q[1];
sx q[1];
rz(0.020843055) q[1];
x q[2];
rz(1.2233569) q[3];
sx q[3];
rz(-1.8104324) q[3];
sx q[3];
rz(0.51260603) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.0032234) q[2];
sx q[2];
rz(-2.1681163) q[2];
sx q[2];
rz(-1.5039911) q[2];
rz(-1.9231632) q[3];
sx q[3];
rz(-2.7259493) q[3];
sx q[3];
rz(0.98201069) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
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
rz(1.0729436) q[0];
sx q[0];
rz(-1.8953841) q[0];
sx q[0];
rz(2.9856227) q[0];
rz(2.7929557) q[1];
sx q[1];
rz(-0.60395423) q[1];
sx q[1];
rz(-1.2329996) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1240163) q[0];
sx q[0];
rz(-1.6917545) q[0];
sx q[0];
rz(1.308754) q[0];
x q[1];
rz(0.92278752) q[2];
sx q[2];
rz(-1.4013259) q[2];
sx q[2];
rz(3.0344935) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.2501331) q[1];
sx q[1];
rz(-0.95055184) q[1];
sx q[1];
rz(2.7930082) q[1];
rz(-pi) q[2];
rz(-0.25147049) q[3];
sx q[3];
rz(-1.5565539) q[3];
sx q[3];
rz(2.3492658) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.6984581) q[2];
sx q[2];
rz(-0.036018697) q[2];
sx q[2];
rz(1.5114463) q[2];
rz(1.1394507) q[3];
sx q[3];
rz(-1.3316493) q[3];
sx q[3];
rz(-2.1512234) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.776942) q[0];
sx q[0];
rz(-2.7237837) q[0];
sx q[0];
rz(1.1530217) q[0];
rz(0.54368377) q[1];
sx q[1];
rz(-1.2666356) q[1];
sx q[1];
rz(2.8797454) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1516929) q[0];
sx q[0];
rz(-1.6542871) q[0];
sx q[0];
rz(-3.1079453) q[0];
rz(-pi) q[1];
rz(-1.0864429) q[2];
sx q[2];
rz(-2.5005955) q[2];
sx q[2];
rz(0.5400368) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.8229554) q[1];
sx q[1];
rz(-1.793211) q[1];
sx q[1];
rz(2.5435124) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.366794) q[3];
sx q[3];
rz(-0.79909924) q[3];
sx q[3];
rz(2.4287756) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5220945) q[2];
sx q[2];
rz(-2.5358989) q[2];
sx q[2];
rz(1.7363133) q[2];
rz(-2.3960579) q[3];
sx q[3];
rz(-1.2657974) q[3];
sx q[3];
rz(2.1982927) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0010506823) q[0];
sx q[0];
rz(-0.11740919) q[0];
sx q[0];
rz(-0.35181272) q[0];
rz(2.70128) q[1];
sx q[1];
rz(-1.3654717) q[1];
sx q[1];
rz(-0.17131677) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85799828) q[0];
sx q[0];
rz(-0.71526113) q[0];
sx q[0];
rz(0.37566988) q[0];
rz(0.75320737) q[2];
sx q[2];
rz(-1.79617) q[2];
sx q[2];
rz(1.6373375) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3371463) q[1];
sx q[1];
rz(-2.4427572) q[1];
sx q[1];
rz(-0.032445907) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2619352) q[3];
sx q[3];
rz(-0.98239693) q[3];
sx q[3];
rz(0.0042875687) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-3.0773086) q[2];
sx q[2];
rz(-1.7694387) q[2];
sx q[2];
rz(-2.1997931) q[2];
rz(1.1549548) q[3];
sx q[3];
rz(-0.20808163) q[3];
sx q[3];
rz(-1.340516) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44523859) q[0];
sx q[0];
rz(-1.1055163) q[0];
sx q[0];
rz(-3.136694) q[0];
rz(-2.694963) q[1];
sx q[1];
rz(-2.547867) q[1];
sx q[1];
rz(-0.68797025) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.60058355) q[0];
sx q[0];
rz(-0.86555153) q[0];
sx q[0];
rz(-1.7974822) q[0];
x q[1];
rz(2.5227658) q[2];
sx q[2];
rz(-1.8558673) q[2];
sx q[2];
rz(3.0322591) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.35536218) q[1];
sx q[1];
rz(-1.3578102) q[1];
sx q[1];
rz(-1.1287685) q[1];
x q[2];
rz(2.3456081) q[3];
sx q[3];
rz(-0.86382691) q[3];
sx q[3];
rz(0.10315264) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.61116162) q[2];
sx q[2];
rz(-0.40803424) q[2];
sx q[2];
rz(2.2116057) q[2];
rz(1.9200578) q[3];
sx q[3];
rz(-2.0285716) q[3];
sx q[3];
rz(2.5482224) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68194836) q[0];
sx q[0];
rz(-1.3197897) q[0];
sx q[0];
rz(0.090106877) q[0];
rz(-1.4247165) q[1];
sx q[1];
rz(-2.5017891) q[1];
sx q[1];
rz(-2.7065014) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2887569) q[0];
sx q[0];
rz(-0.52271862) q[0];
sx q[0];
rz(-2.8466892) q[0];
rz(2.4736011) q[2];
sx q[2];
rz(-1.4282303) q[2];
sx q[2];
rz(0.1145471) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.7104946) q[1];
sx q[1];
rz(-2.5113547) q[1];
sx q[1];
rz(1.3517387) q[1];
rz(-pi) q[2];
rz(3.1396542) q[3];
sx q[3];
rz(-0.42436212) q[3];
sx q[3];
rz(0.75943179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.6626176) q[2];
sx q[2];
rz(-2.0765897) q[2];
sx q[2];
rz(2.5153861) q[2];
rz(2.9446757) q[3];
sx q[3];
rz(-0.99072376) q[3];
sx q[3];
rz(-1.3885434) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5597252) q[0];
sx q[0];
rz(-1.6904866) q[0];
sx q[0];
rz(-3.0873121) q[0];
rz(0.42298969) q[1];
sx q[1];
rz(-1.3754247) q[1];
sx q[1];
rz(-1.3351701) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2535219) q[0];
sx q[0];
rz(-1.0253064) q[0];
sx q[0];
rz(-2.2175199) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0112052) q[2];
sx q[2];
rz(-1.2647243) q[2];
sx q[2];
rz(-2.0843506) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.21653342) q[1];
sx q[1];
rz(-2.0833942) q[1];
sx q[1];
rz(-2.53043) q[1];
x q[2];
rz(1.1109675) q[3];
sx q[3];
rz(-0.15367344) q[3];
sx q[3];
rz(2.3914162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.52428952) q[2];
sx q[2];
rz(-1.6733988) q[2];
sx q[2];
rz(-0.35150251) q[2];
rz(1.7678123) q[3];
sx q[3];
rz(-0.51353729) q[3];
sx q[3];
rz(2.8721749) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7518625) q[0];
sx q[0];
rz(-0.26079145) q[0];
sx q[0];
rz(-0.75916284) q[0];
rz(1.0836541) q[1];
sx q[1];
rz(-1.5307129) q[1];
sx q[1];
rz(-1.0677451) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5006951) q[0];
sx q[0];
rz(-1.0984549) q[0];
sx q[0];
rz(-0.91820902) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.75394404) q[2];
sx q[2];
rz(-0.39421668) q[2];
sx q[2];
rz(-2.8682402) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.8651004) q[1];
sx q[1];
rz(-2.4363764) q[1];
sx q[1];
rz(0.92030763) q[1];
rz(-pi) q[2];
rz(1.9633406) q[3];
sx q[3];
rz(-2.4053229) q[3];
sx q[3];
rz(-0.53680778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.4296253) q[2];
sx q[2];
rz(-1.2099313) q[2];
sx q[2];
rz(0.21027002) q[2];
rz(-0.78768864) q[3];
sx q[3];
rz(-1.382788) q[3];
sx q[3];
rz(0.53019607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6982211) q[0];
sx q[0];
rz(-2.5826695) q[0];
sx q[0];
rz(-1.2179751) q[0];
rz(1.632985) q[1];
sx q[1];
rz(-1.8764381) q[1];
sx q[1];
rz(1.1988342) q[1];
rz(2.2589113) q[2];
sx q[2];
rz(-1.8919049) q[2];
sx q[2];
rz(1.6406825) q[2];
rz(1.2959215) q[3];
sx q[3];
rz(-1.4209005) q[3];
sx q[3];
rz(3.0358547) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
