OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.6141619) q[0];
sx q[0];
rz(-0.80978137) q[0];
sx q[0];
rz(-0.53139395) q[0];
rz(0.2375138) q[1];
sx q[1];
rz(1.3637435) q[1];
sx q[1];
rz(10.627828) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0396084) q[0];
sx q[0];
rz(-1.7829722) q[0];
sx q[0];
rz(2.9456445) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8868944) q[2];
sx q[2];
rz(-1.7809976) q[2];
sx q[2];
rz(2.3496698) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.4431475) q[1];
sx q[1];
rz(-1.5138211) q[1];
sx q[1];
rz(0.2351825) q[1];
rz(-0.12227998) q[3];
sx q[3];
rz(-2.8465726) q[3];
sx q[3];
rz(-0.89000851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.27665859) q[2];
sx q[2];
rz(-1.1316789) q[2];
sx q[2];
rz(-0.19908389) q[2];
rz(-1.52786) q[3];
sx q[3];
rz(-0.49694967) q[3];
sx q[3];
rz(-0.73408192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1428225) q[0];
sx q[0];
rz(-0.12717371) q[0];
sx q[0];
rz(2.3221827) q[0];
rz(-0.2858513) q[1];
sx q[1];
rz(-2.2276623) q[1];
sx q[1];
rz(-1.2664638) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3931912) q[0];
sx q[0];
rz(-1.1921492) q[0];
sx q[0];
rz(-2.4875531) q[0];
rz(-2.9026047) q[2];
sx q[2];
rz(-1.1766489) q[2];
sx q[2];
rz(0.1220526) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.30217583) q[1];
sx q[1];
rz(-1.2100879) q[1];
sx q[1];
rz(-1.9990666) q[1];
rz(-2.0224138) q[3];
sx q[3];
rz(-1.041881) q[3];
sx q[3];
rz(0.12714735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.1566029) q[2];
sx q[2];
rz(-1.2133602) q[2];
sx q[2];
rz(-1.161969) q[2];
rz(0.087163838) q[3];
sx q[3];
rz(-1.6379084) q[3];
sx q[3];
rz(-0.073908977) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53482985) q[0];
sx q[0];
rz(-2.02089) q[0];
sx q[0];
rz(2.8379922) q[0];
rz(1.7595694) q[1];
sx q[1];
rz(-1.2761812) q[1];
sx q[1];
rz(0.64750013) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.835639) q[0];
sx q[0];
rz(-1.5257972) q[0];
sx q[0];
rz(2.5070028) q[0];
x q[1];
rz(2.5325003) q[2];
sx q[2];
rz(-2.3272772) q[2];
sx q[2];
rz(0.44037214) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.68413098) q[1];
sx q[1];
rz(-2.0254454) q[1];
sx q[1];
rz(0.40747868) q[1];
rz(-pi) q[2];
rz(0.89214274) q[3];
sx q[3];
rz(-2.3193079) q[3];
sx q[3];
rz(3.126614) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.8422164) q[2];
sx q[2];
rz(-2.3589578) q[2];
sx q[2];
rz(1.5608609) q[2];
rz(-2.1598699) q[3];
sx q[3];
rz(-1.2795307) q[3];
sx q[3];
rz(2.4981807) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9075539) q[0];
sx q[0];
rz(-2.3038395) q[0];
sx q[0];
rz(-2.2696944) q[0];
rz(0.38726989) q[1];
sx q[1];
rz(-0.64278066) q[1];
sx q[1];
rz(-2.8505468) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0204578) q[0];
sx q[0];
rz(-0.99636787) q[0];
sx q[0];
rz(-2.8681884) q[0];
rz(-pi) q[1];
rz(-1.8918858) q[2];
sx q[2];
rz(-1.2920657) q[2];
sx q[2];
rz(2.2433787) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.0709666) q[1];
sx q[1];
rz(-0.91849209) q[1];
sx q[1];
rz(0.91826203) q[1];
rz(0.16104161) q[3];
sx q[3];
rz(-2.1412009) q[3];
sx q[3];
rz(2.963484) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.7724093) q[2];
sx q[2];
rz(-2.3466551) q[2];
sx q[2];
rz(2.5975361) q[2];
rz(-0.50928515) q[3];
sx q[3];
rz(-2.0777168) q[3];
sx q[3];
rz(0.86597401) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0020224) q[0];
sx q[0];
rz(-2.1152571) q[0];
sx q[0];
rz(2.496526) q[0];
rz(2.5158665) q[1];
sx q[1];
rz(-1.2236243) q[1];
sx q[1];
rz(-2.5114139) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7932927) q[0];
sx q[0];
rz(-0.59069809) q[0];
sx q[0];
rz(-0.22038711) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1516045) q[2];
sx q[2];
rz(-2.345511) q[2];
sx q[2];
rz(0.58005709) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.6781792) q[1];
sx q[1];
rz(-2.5555875) q[1];
sx q[1];
rz(1.9834118) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.583902) q[3];
sx q[3];
rz(-1.7002749) q[3];
sx q[3];
rz(1.492471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.75382918) q[2];
sx q[2];
rz(-1.3239653) q[2];
sx q[2];
rz(1.7774263) q[2];
rz(1.8917313) q[3];
sx q[3];
rz(-1.9259689) q[3];
sx q[3];
rz(2.2201339) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0867778) q[0];
sx q[0];
rz(-2.5578816) q[0];
sx q[0];
rz(-2.4247647) q[0];
rz(-2.7038799) q[1];
sx q[1];
rz(-2.4611459) q[1];
sx q[1];
rz(0.81370083) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1610634) q[0];
sx q[0];
rz(-1.6335532) q[0];
sx q[0];
rz(-1.3570696) q[0];
rz(-pi) q[1];
rz(0.99102334) q[2];
sx q[2];
rz(-0.93587854) q[2];
sx q[2];
rz(2.6099176) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.70255792) q[1];
sx q[1];
rz(-0.82779373) q[1];
sx q[1];
rz(2.031416) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2171451) q[3];
sx q[3];
rz(-1.2815223) q[3];
sx q[3];
rz(-1.7464964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.2563236) q[2];
sx q[2];
rz(-2.7041114) q[2];
sx q[2];
rz(-1.1876594) q[2];
rz(0.93196431) q[3];
sx q[3];
rz(-2.0075802) q[3];
sx q[3];
rz(0.37117547) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9880144) q[0];
sx q[0];
rz(-1.0289285) q[0];
sx q[0];
rz(2.4555092) q[0];
rz(1.6185121) q[1];
sx q[1];
rz(-2.1673817) q[1];
sx q[1];
rz(0.66326052) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2880741) q[0];
sx q[0];
rz(-0.92068866) q[0];
sx q[0];
rz(3.1115565) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3528321) q[2];
sx q[2];
rz(-0.94467615) q[2];
sx q[2];
rz(1.2499714) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.82603964) q[1];
sx q[1];
rz(-1.3693046) q[1];
sx q[1];
rz(-2.0622581) q[1];
rz(2.0459941) q[3];
sx q[3];
rz(-1.0605269) q[3];
sx q[3];
rz(-2.1197532) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.66701165) q[2];
sx q[2];
rz(-2.3040999) q[2];
sx q[2];
rz(-0.015080301) q[2];
rz(1.0117426) q[3];
sx q[3];
rz(-1.116131) q[3];
sx q[3];
rz(0.57141465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.052977234) q[0];
sx q[0];
rz(-0.73260728) q[0];
sx q[0];
rz(0.10738871) q[0];
rz(-0.30474162) q[1];
sx q[1];
rz(-1.823103) q[1];
sx q[1];
rz(-2.6838141) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7664238) q[0];
sx q[0];
rz(-0.94978118) q[0];
sx q[0];
rz(0.39874886) q[0];
rz(-pi) q[1];
rz(-0.56577487) q[2];
sx q[2];
rz(-0.83116311) q[2];
sx q[2];
rz(0.6643675) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.0486054) q[1];
sx q[1];
rz(-0.43097207) q[1];
sx q[1];
rz(-0.98577568) q[1];
rz(-pi) q[2];
rz(2.5910208) q[3];
sx q[3];
rz(-0.62567657) q[3];
sx q[3];
rz(1.613137) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7765939) q[2];
sx q[2];
rz(-1.685131) q[2];
sx q[2];
rz(1.7101074) q[2];
rz(1.4252023) q[3];
sx q[3];
rz(-2.5148354) q[3];
sx q[3];
rz(1.2861929) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15437056) q[0];
sx q[0];
rz(-0.4796589) q[0];
sx q[0];
rz(-0.95348683) q[0];
rz(1.8158688) q[1];
sx q[1];
rz(-0.49294254) q[1];
sx q[1];
rz(-2.5568533) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.23823243) q[0];
sx q[0];
rz(-2.0582709) q[0];
sx q[0];
rz(-1.4438629) q[0];
x q[1];
rz(-2.598001) q[2];
sx q[2];
rz(-1.2541176) q[2];
sx q[2];
rz(1.2515595) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.1251724) q[1];
sx q[1];
rz(-0.44508815) q[1];
sx q[1];
rz(-0.16146407) q[1];
rz(1.3034079) q[3];
sx q[3];
rz(-1.7087473) q[3];
sx q[3];
rz(1.4064058) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.3328302) q[2];
sx q[2];
rz(-0.82020438) q[2];
sx q[2];
rz(2.8653223) q[2];
rz(2.590498) q[3];
sx q[3];
rz(-1.7421744) q[3];
sx q[3];
rz(2.7579257) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0722512) q[0];
sx q[0];
rz(-2.6273917) q[0];
sx q[0];
rz(-0.65548354) q[0];
rz(1.1570702) q[1];
sx q[1];
rz(-1.7600704) q[1];
sx q[1];
rz(1.5225333) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5009506) q[0];
sx q[0];
rz(-2.2146533) q[0];
sx q[0];
rz(1.9615016) q[0];
x q[1];
rz(-2.7340639) q[2];
sx q[2];
rz(-0.73020836) q[2];
sx q[2];
rz(-0.66116316) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.9177502) q[1];
sx q[1];
rz(-2.192944) q[1];
sx q[1];
rz(-1.7734581) q[1];
rz(-2.7195815) q[3];
sx q[3];
rz(-1.9415628) q[3];
sx q[3];
rz(-1.5164204) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.15554252) q[2];
sx q[2];
rz(-2.7337998) q[2];
sx q[2];
rz(-0.84214169) q[2];
rz(-1.5367674) q[3];
sx q[3];
rz(-1.3804881) q[3];
sx q[3];
rz(-3.1150637) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90606541) q[0];
sx q[0];
rz(-1.9530095) q[0];
sx q[0];
rz(-0.45146913) q[0];
rz(-2.4189667) q[1];
sx q[1];
rz(-0.87651785) q[1];
sx q[1];
rz(-1.6095907) q[1];
rz(-0.026635219) q[2];
sx q[2];
rz(-1.7799108) q[2];
sx q[2];
rz(-1.5540661) q[2];
rz(1.6668672) q[3];
sx q[3];
rz(-2.4484652) q[3];
sx q[3];
rz(-0.0948003) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
