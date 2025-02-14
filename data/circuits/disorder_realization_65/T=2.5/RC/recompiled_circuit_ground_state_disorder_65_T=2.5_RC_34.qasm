OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.559691) q[0];
sx q[0];
rz(-0.087495916) q[0];
sx q[0];
rz(-3.0425332) q[0];
rz(-0.29095185) q[1];
sx q[1];
rz(4.8922242) q[1];
sx q[1];
rz(7.7982245) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6884424) q[0];
sx q[0];
rz(-0.52061284) q[0];
sx q[0];
rz(-2.0474035) q[0];
rz(0.85364437) q[2];
sx q[2];
rz(-0.63736299) q[2];
sx q[2];
rz(1.2946707) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.6219306) q[1];
sx q[1];
rz(-2.3708713) q[1];
sx q[1];
rz(1.3407441) q[1];
rz(-pi) q[2];
x q[2];
rz(0.73554773) q[3];
sx q[3];
rz(-0.41829073) q[3];
sx q[3];
rz(-2.4052704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0809975) q[2];
sx q[2];
rz(-0.011757714) q[2];
sx q[2];
rz(-0.11059977) q[2];
rz(-3.1317173) q[3];
sx q[3];
rz(-0.014009352) q[3];
sx q[3];
rz(-2.2601155) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1572384) q[0];
sx q[0];
rz(-0.088033661) q[0];
sx q[0];
rz(1.3798168) q[0];
rz(3.1294322) q[1];
sx q[1];
rz(-2.1575243) q[1];
sx q[1];
rz(-1.5252569) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75689689) q[0];
sx q[0];
rz(-1.257557) q[0];
sx q[0];
rz(-1.6301073) q[0];
rz(-pi) q[1];
x q[1];
rz(0.42945822) q[2];
sx q[2];
rz(-2.6864687) q[2];
sx q[2];
rz(-1.5433102) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-3.1327019) q[1];
sx q[1];
rz(-1.3976919) q[1];
sx q[1];
rz(-1.6457896) q[1];
x q[2];
rz(-0.44604519) q[3];
sx q[3];
rz(-1.6823263) q[3];
sx q[3];
rz(-2.1457637) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.81185594) q[2];
sx q[2];
rz(-0.024710329) q[2];
sx q[2];
rz(-0.032026637) q[2];
rz(0.61102593) q[3];
sx q[3];
rz(-1.150584) q[3];
sx q[3];
rz(1.7648296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5585612) q[0];
sx q[0];
rz(-0.91201454) q[0];
sx q[0];
rz(1.6770021) q[0];
rz(-1.6327845) q[1];
sx q[1];
rz(-1.8784411) q[1];
sx q[1];
rz(0.060922932) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7446242) q[0];
sx q[0];
rz(-1.3875657) q[0];
sx q[0];
rz(0.38559648) q[0];
rz(-pi) q[1];
rz(1.6834358) q[2];
sx q[2];
rz(-1.4667744) q[2];
sx q[2];
rz(-2.011781) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(3.0197501) q[1];
sx q[1];
rz(-1.3977244) q[1];
sx q[1];
rz(-0.16000859) q[1];
rz(-pi) q[2];
rz(-2.4931566) q[3];
sx q[3];
rz(-0.47659527) q[3];
sx q[3];
rz(2.4364382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.9691201) q[2];
sx q[2];
rz(-3.1388333) q[2];
sx q[2];
rz(0.74392444) q[2];
rz(-0.38532358) q[3];
sx q[3];
rz(-3.139956) q[3];
sx q[3];
rz(2.2374432) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1291523) q[0];
sx q[0];
rz(-3.0991982) q[0];
sx q[0];
rz(2.2989035) q[0];
rz(2.1878302) q[1];
sx q[1];
rz(-1.5019491) q[1];
sx q[1];
rz(-1.5488254) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7208389) q[0];
sx q[0];
rz(-0.11891236) q[0];
sx q[0];
rz(-0.17158385) q[0];
rz(-pi) q[1];
rz(1.2090525) q[2];
sx q[2];
rz(-0.42297034) q[2];
sx q[2];
rz(1.376938) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.4348244) q[1];
sx q[1];
rz(-2.206496) q[1];
sx q[1];
rz(3.0830361) q[1];
rz(-pi) q[2];
x q[2];
rz(0.10080849) q[3];
sx q[3];
rz(-0.78323302) q[3];
sx q[3];
rz(0.99683652) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-3.039975) q[2];
sx q[2];
rz(-3.1246694) q[2];
sx q[2];
rz(-0.42538154) q[2];
rz(1.7915122) q[3];
sx q[3];
rz(-1.5964419) q[3];
sx q[3];
rz(2.8369956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4199583) q[0];
sx q[0];
rz(-3.1379134) q[0];
sx q[0];
rz(-0.71596181) q[0];
rz(1.51659) q[1];
sx q[1];
rz(-0.45418987) q[1];
sx q[1];
rz(2.9862278) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6139406) q[0];
sx q[0];
rz(-0.10738534) q[0];
sx q[0];
rz(-1.7201109) q[0];
rz(-pi) q[1];
rz(2.7648126) q[2];
sx q[2];
rz(-1.6544806) q[2];
sx q[2];
rz(2.7629536) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.027861431) q[1];
sx q[1];
rz(-1.3044954) q[1];
sx q[1];
rz(1.6990425) q[1];
x q[2];
rz(0.89837667) q[3];
sx q[3];
rz(-1.9168609) q[3];
sx q[3];
rz(-2.4806674) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.4510437) q[2];
sx q[2];
rz(-2.0205708) q[2];
sx q[2];
rz(1.5345908) q[2];
rz(2.0524041) q[3];
sx q[3];
rz(-0.023622731) q[3];
sx q[3];
rz(-1.7781809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4264193) q[0];
sx q[0];
rz(-0.024113163) q[0];
sx q[0];
rz(0.56060785) q[0];
rz(2.2757065) q[1];
sx q[1];
rz(-0.0090323369) q[1];
sx q[1];
rz(2.3057356) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27640057) q[0];
sx q[0];
rz(-1.5888583) q[0];
sx q[0];
rz(-2.8457321) q[0];
rz(-pi) q[1];
rz(-1.1634665) q[2];
sx q[2];
rz(-1.6025123) q[2];
sx q[2];
rz(-2.9518366) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.8053956) q[1];
sx q[1];
rz(-1.312689) q[1];
sx q[1];
rz(1.6805524) q[1];
rz(1.3243598) q[3];
sx q[3];
rz(-1.8772238) q[3];
sx q[3];
rz(0.6014733) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.5777638) q[2];
sx q[2];
rz(-1.8356297) q[2];
sx q[2];
rz(-1.5527661) q[2];
rz(1.0519489) q[3];
sx q[3];
rz(-0.51783872) q[3];
sx q[3];
rz(0.53798211) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2688667) q[0];
sx q[0];
rz(-2.5476542) q[0];
sx q[0];
rz(-1.8178222) q[0];
rz(2.2920604) q[1];
sx q[1];
rz(-3.140026) q[1];
sx q[1];
rz(0.82226396) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0104843) q[0];
sx q[0];
rz(-2.9379651) q[0];
sx q[0];
rz(-2.0631316) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9862764) q[2];
sx q[2];
rz(-0.34201038) q[2];
sx q[2];
rz(-1.6928455) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.7436372) q[1];
sx q[1];
rz(-1.5642867) q[1];
sx q[1];
rz(1.5012451) q[1];
rz(-pi) q[2];
rz(2.5093171) q[3];
sx q[3];
rz(-0.32514363) q[3];
sx q[3];
rz(2.7141311) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.5620586) q[2];
sx q[2];
rz(-1.2349393) q[2];
sx q[2];
rz(2.298992) q[2];
rz(-0.54251999) q[3];
sx q[3];
rz(-0.018535651) q[3];
sx q[3];
rz(0.59244573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0845959) q[0];
sx q[0];
rz(-1.5472941) q[0];
sx q[0];
rz(-2.090276) q[0];
rz(-1.7683138) q[1];
sx q[1];
rz(-0.028742464) q[1];
sx q[1];
rz(1.2398667) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43869685) q[0];
sx q[0];
rz(-1.6099766) q[0];
sx q[0];
rz(0.24287276) q[0];
x q[1];
rz(-1.1781111) q[2];
sx q[2];
rz(-1.6047021) q[2];
sx q[2];
rz(0.086080649) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(3.0770331) q[1];
sx q[1];
rz(-1.7747595) q[1];
sx q[1];
rz(2.166421) q[1];
x q[2];
rz(-1.3040335) q[3];
sx q[3];
rz(-1.3893616) q[3];
sx q[3];
rz(0.72600466) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.54590589) q[2];
sx q[2];
rz(-0.02820708) q[2];
sx q[2];
rz(-2.9941881) q[2];
rz(1.5386511) q[3];
sx q[3];
rz(-1.4378865) q[3];
sx q[3];
rz(-2.9586207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.962917) q[0];
sx q[0];
rz(-2.8067532) q[0];
sx q[0];
rz(-1.8033002) q[0];
rz(-1.6637404) q[1];
sx q[1];
rz(-2.0070952) q[1];
sx q[1];
rz(-2.9683364) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4497183) q[0];
sx q[0];
rz(-2.1341748) q[0];
sx q[0];
rz(2.0823576) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1384754) q[2];
sx q[2];
rz(-2.4698493) q[2];
sx q[2];
rz(0.22638182) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4819255) q[1];
sx q[1];
rz(-0.29259455) q[1];
sx q[1];
rz(-0.56168075) q[1];
rz(-pi) q[2];
rz(-1.5842755) q[3];
sx q[3];
rz(-1.5649678) q[3];
sx q[3];
rz(-1.7268501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.0080537) q[2];
sx q[2];
rz(-3.1296788) q[2];
sx q[2];
rz(-2.1214205) q[2];
rz(-0.080848761) q[3];
sx q[3];
rz(-0.0011708502) q[3];
sx q[3];
rz(1.1297869) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5197649) q[0];
sx q[0];
rz(-2.8563359) q[0];
sx q[0];
rz(1.4813625) q[0];
rz(0.18352428) q[1];
sx q[1];
rz(-0.32569519) q[1];
sx q[1];
rz(0.087800177) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7526898) q[0];
sx q[0];
rz(-1.6713033) q[0];
sx q[0];
rz(-0.44403063) q[0];
rz(-pi) q[1];
rz(-1.4884639) q[2];
sx q[2];
rz(-1.265324) q[2];
sx q[2];
rz(2.9412172) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.99276772) q[1];
sx q[1];
rz(-1.449905) q[1];
sx q[1];
rz(2.8962027) q[1];
rz(2.4724768) q[3];
sx q[3];
rz(-1.8103616) q[3];
sx q[3];
rz(-1.4544044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.2878908) q[2];
sx q[2];
rz(-0.012306865) q[2];
sx q[2];
rz(2.2575209) q[2];
rz(-2.4331802) q[3];
sx q[3];
rz(-3.1413779) q[3];
sx q[3];
rz(-0.29402548) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52083279) q[0];
sx q[0];
rz(-1.5694869) q[0];
sx q[0];
rz(-1.6318305) q[0];
rz(-3.1132501) q[1];
sx q[1];
rz(-0.62348532) q[1];
sx q[1];
rz(-3.0709406) q[1];
rz(-2.6940343) q[2];
sx q[2];
rz(-1.9050514) q[2];
sx q[2];
rz(-1.5290268) q[2];
rz(0.36793637) q[3];
sx q[3];
rz(-0.76303861) q[3];
sx q[3];
rz(-1.7581024) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
