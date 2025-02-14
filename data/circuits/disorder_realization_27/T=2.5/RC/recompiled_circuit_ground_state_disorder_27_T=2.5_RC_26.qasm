OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.55750027) q[0];
sx q[0];
rz(-3.1182365) q[0];
sx q[0];
rz(-0.93552247) q[0];
rz(-1.6159396) q[1];
sx q[1];
rz(-1.5204117) q[1];
sx q[1];
rz(2.867155) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4545499) q[0];
sx q[0];
rz(-1.4996212) q[0];
sx q[0];
rz(0.068333621) q[0];
rz(-pi) q[1];
rz(-0.69752946) q[2];
sx q[2];
rz(-2.4781057) q[2];
sx q[2];
rz(1.0349719) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1573616) q[1];
sx q[1];
rz(-1.5604291) q[1];
sx q[1];
rz(-1.5342516) q[1];
x q[2];
rz(1.4001582) q[3];
sx q[3];
rz(-1.4981761) q[3];
sx q[3];
rz(-0.69334465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.9258257) q[2];
sx q[2];
rz(-3.1302858) q[2];
sx q[2];
rz(2.0634148) q[2];
rz(0.82873851) q[3];
sx q[3];
rz(-1.6308035) q[3];
sx q[3];
rz(0.80882788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.082212903) q[0];
sx q[0];
rz(-1.2780715) q[0];
sx q[0];
rz(1.3905806) q[0];
rz(1.4311721) q[1];
sx q[1];
rz(-0.0043967604) q[1];
sx q[1];
rz(1.4336525) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19884685) q[0];
sx q[0];
rz(-2.2042243) q[0];
sx q[0];
rz(-0.1300227) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.016489224) q[2];
sx q[2];
rz(-1.1270583) q[2];
sx q[2];
rz(-1.5373785) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.2428618) q[1];
sx q[1];
rz(-3.1325097) q[1];
sx q[1];
rz(1.6159978) q[1];
rz(-2.3539813) q[3];
sx q[3];
rz(-1.5244841) q[3];
sx q[3];
rz(-2.5758343) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7626875) q[2];
sx q[2];
rz(-1.5451558) q[2];
sx q[2];
rz(1.572466) q[2];
rz(-2.3804741) q[3];
sx q[3];
rz(-3.0877536) q[3];
sx q[3];
rz(1.2083453) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.222027) q[0];
sx q[0];
rz(-2.5313105) q[0];
sx q[0];
rz(-0.56184226) q[0];
rz(1.5737083) q[1];
sx q[1];
rz(-1.578873) q[1];
sx q[1];
rz(0.017339658) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2527232) q[0];
sx q[0];
rz(-2.4742352) q[0];
sx q[0];
rz(-0.86020893) q[0];
rz(-pi) q[1];
rz(2.1465535) q[2];
sx q[2];
rz(-0.77313609) q[2];
sx q[2];
rz(0.38627689) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.8120809) q[1];
sx q[1];
rz(-1.2083665) q[1];
sx q[1];
rz(-3.1245438) q[1];
rz(-pi) q[2];
rz(3.1005834) q[3];
sx q[3];
rz(-2.6309391) q[3];
sx q[3];
rz(1.2488431) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.5138381) q[2];
sx q[2];
rz(-1.4189812) q[2];
sx q[2];
rz(-0.55297744) q[2];
rz(-1.2051469) q[3];
sx q[3];
rz(-1.5744934) q[3];
sx q[3];
rz(-1.5945826) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39128458) q[0];
sx q[0];
rz(-2.1109844) q[0];
sx q[0];
rz(-1.3543825) q[0];
rz(1.7693819) q[1];
sx q[1];
rz(-0.0028227614) q[1];
sx q[1];
rz(1.7599958) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.87028044) q[0];
sx q[0];
rz(-0.63374937) q[0];
sx q[0];
rz(-0.84931727) q[0];
rz(-pi) q[1];
rz(-2.4468378) q[2];
sx q[2];
rz(-3.1379897) q[2];
sx q[2];
rz(-2.3191593) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.957037) q[1];
sx q[1];
rz(-0.77691764) q[1];
sx q[1];
rz(-2.0691815) q[1];
rz(0.28691157) q[3];
sx q[3];
rz(-2.3182456) q[3];
sx q[3];
rz(0.2940184) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.0733205) q[2];
sx q[2];
rz(-0.019651532) q[2];
sx q[2];
rz(1.9015296) q[2];
rz(2.9114919) q[3];
sx q[3];
rz(-3.1374044) q[3];
sx q[3];
rz(0.41964644) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0536026) q[0];
sx q[0];
rz(-0.78730655) q[0];
sx q[0];
rz(-1.8863652) q[0];
rz(-0.0093731006) q[1];
sx q[1];
rz(-1.3690989) q[1];
sx q[1];
rz(3.110041) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9242212) q[0];
sx q[0];
rz(-1.517202) q[0];
sx q[0];
rz(1.0384667) q[0];
x q[1];
rz(-1.8589694) q[2];
sx q[2];
rz(-2.8409578) q[2];
sx q[2];
rz(0.86821454) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.0134558) q[1];
sx q[1];
rz(-1.0425048) q[1];
sx q[1];
rz(2.8847413) q[1];
x q[2];
rz(-2.3208951) q[3];
sx q[3];
rz(-2.3682171) q[3];
sx q[3];
rz(3.000963) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.3343398) q[2];
sx q[2];
rz(-0.0060609239) q[2];
sx q[2];
rz(-0.80017153) q[2];
rz(2.3470894) q[3];
sx q[3];
rz(-0.032363351) q[3];
sx q[3];
rz(2.2728424) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0873347) q[0];
sx q[0];
rz(-0.15905173) q[0];
sx q[0];
rz(1.4999088) q[0];
rz(2.9686887) q[1];
sx q[1];
rz(-0.042363107) q[1];
sx q[1];
rz(-0.049887966) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5373851) q[0];
sx q[0];
rz(-2.4695463) q[0];
sx q[0];
rz(1.1623357) q[0];
rz(-pi) q[1];
rz(-0.11325963) q[2];
sx q[2];
rz(-2.335603) q[2];
sx q[2];
rz(1.9085371) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.7916471) q[1];
sx q[1];
rz(-2.269372) q[1];
sx q[1];
rz(2.8142745) q[1];
rz(-0.84953725) q[3];
sx q[3];
rz(-1.153933) q[3];
sx q[3];
rz(0.45826926) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.738203) q[2];
sx q[2];
rz(-3.093779) q[2];
sx q[2];
rz(-1.2460463) q[2];
rz(-1.3561148) q[3];
sx q[3];
rz(-0.035364371) q[3];
sx q[3];
rz(1.416392) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7127011) q[0];
sx q[0];
rz(-2.3172947) q[0];
sx q[0];
rz(-1.7190546) q[0];
rz(1.9046344) q[1];
sx q[1];
rz(-0.037038602) q[1];
sx q[1];
rz(-0.2027771) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0109459) q[0];
sx q[0];
rz(-1.2460862) q[0];
sx q[0];
rz(-0.79239158) q[0];
x q[1];
rz(-2.3788664) q[2];
sx q[2];
rz(-0.87050754) q[2];
sx q[2];
rz(1.8000079) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.4331268) q[1];
sx q[1];
rz(-0.34516343) q[1];
sx q[1];
rz(0.16394798) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0596541) q[3];
sx q[3];
rz(-1.7082001) q[3];
sx q[3];
rz(-2.5554772) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.612959) q[2];
sx q[2];
rz(-0.10047675) q[2];
sx q[2];
rz(-0.67451492) q[2];
rz(1.3450735) q[3];
sx q[3];
rz(-2.9967872) q[3];
sx q[3];
rz(-1.8176746) q[3];
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
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26789185) q[0];
sx q[0];
rz(-0.75650263) q[0];
sx q[0];
rz(0.85195136) q[0];
rz(-0.1918699) q[1];
sx q[1];
rz(-3.1286897) q[1];
sx q[1];
rz(2.875944) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5776731) q[0];
sx q[0];
rz(-1.1611337) q[0];
sx q[0];
rz(1.3807317) q[0];
x q[1];
rz(2.526409) q[2];
sx q[2];
rz(-1.8447723) q[2];
sx q[2];
rz(0.09466234) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.3623464) q[1];
sx q[1];
rz(-1.4873449) q[1];
sx q[1];
rz(-1.6298184) q[1];
rz(-0.55691584) q[3];
sx q[3];
rz(-1.89291) q[3];
sx q[3];
rz(0.39469257) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.3702281) q[2];
sx q[2];
rz(-3.0493272) q[2];
sx q[2];
rz(2.6293758) q[2];
rz(0.15906119) q[3];
sx q[3];
rz(-3.1064807) q[3];
sx q[3];
rz(-1.4353282) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2559741) q[0];
sx q[0];
rz(-1.4775448) q[0];
sx q[0];
rz(-1.0783827) q[0];
rz(-1.6400317) q[1];
sx q[1];
rz(-2.9402132) q[1];
sx q[1];
rz(1.5837502) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5756861) q[0];
sx q[0];
rz(-2.4543336) q[0];
sx q[0];
rz(-3.1172196) q[0];
rz(1.0859162) q[2];
sx q[2];
rz(-0.74046248) q[2];
sx q[2];
rz(-1.6892576) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.4068102) q[1];
sx q[1];
rz(-1.5487707) q[1];
sx q[1];
rz(-1.5741482) q[1];
rz(-2.7872378) q[3];
sx q[3];
rz(-2.1956596) q[3];
sx q[3];
rz(-2.5884678) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.25643361) q[2];
sx q[2];
rz(-3.1270471) q[2];
sx q[2];
rz(-0.069615901) q[2];
rz(-3.0749248) q[3];
sx q[3];
rz(-1.0144517) q[3];
sx q[3];
rz(0.67924172) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2081864) q[0];
sx q[0];
rz(-1.2766159) q[0];
sx q[0];
rz(0.25190121) q[0];
rz(-1.6690286) q[1];
sx q[1];
rz(-2.9258969) q[1];
sx q[1];
rz(-3.0676945) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4842) q[0];
sx q[0];
rz(-1.2061242) q[0];
sx q[0];
rz(3.0658989) q[0];
rz(-0.017357512) q[2];
sx q[2];
rz(-1.6021172) q[2];
sx q[2];
rz(0.90085122) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-3.0987948) q[1];
sx q[1];
rz(-0.86314161) q[1];
sx q[1];
rz(1.2513729) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4611887) q[3];
sx q[3];
rz(-0.26016737) q[3];
sx q[3];
rz(-0.55458595) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.4613688) q[2];
sx q[2];
rz(-0.0071439925) q[2];
sx q[2];
rz(0.76772493) q[2];
rz(-1.7420306) q[3];
sx q[3];
rz(-0.00082409516) q[3];
sx q[3];
rz(2.6373533) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5117699) q[0];
sx q[0];
rz(-0.98291021) q[0];
sx q[0];
rz(1.7194189) q[0];
rz(-3.1172251) q[1];
sx q[1];
rz(-0.15932803) q[1];
sx q[1];
rz(-2.9111964) q[1];
rz(-0.38300285) q[2];
sx q[2];
rz(-2.2097864) q[2];
sx q[2];
rz(-1.0351444) q[2];
rz(-1.6679933) q[3];
sx q[3];
rz(-1.9862277) q[3];
sx q[3];
rz(1.338892) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
