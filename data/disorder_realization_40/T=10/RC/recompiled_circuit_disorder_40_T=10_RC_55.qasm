OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.6089132) q[0];
sx q[0];
rz(-0.37663868) q[0];
sx q[0];
rz(-3.0298046) q[0];
rz(-1.4594266) q[1];
sx q[1];
rz(-1.6571801) q[1];
sx q[1];
rz(0.15375528) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35225866) q[0];
sx q[0];
rz(-0.52868045) q[0];
sx q[0];
rz(-2.5369011) q[0];
rz(-pi) q[1];
rz(-1.2674238) q[2];
sx q[2];
rz(-0.25250013) q[2];
sx q[2];
rz(-1.4361824) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.44126025) q[1];
sx q[1];
rz(-1.889519) q[1];
sx q[1];
rz(3.1096427) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1313685) q[3];
sx q[3];
rz(-1.8779552) q[3];
sx q[3];
rz(2.9205521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.6979606) q[2];
sx q[2];
rz(-1.4322832) q[2];
sx q[2];
rz(-1.4367746) q[2];
rz(-0.73389655) q[3];
sx q[3];
rz(-1.5489483) q[3];
sx q[3];
rz(2.6255887) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.88965082) q[0];
sx q[0];
rz(-1.2263068) q[0];
sx q[0];
rz(-0.92457986) q[0];
rz(0.997116) q[1];
sx q[1];
rz(-2.6328502) q[1];
sx q[1];
rz(1.8181713) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2307777) q[0];
sx q[0];
rz(-1.9998904) q[0];
sx q[0];
rz(-0.50165117) q[0];
rz(-pi) q[1];
rz(-0.1214059) q[2];
sx q[2];
rz(-2.7896023) q[2];
sx q[2];
rz(-2.6999203) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0022701) q[1];
sx q[1];
rz(-1.8404669) q[1];
sx q[1];
rz(-0.76489277) q[1];
x q[2];
rz(2.1341392) q[3];
sx q[3];
rz(-2.4881119) q[3];
sx q[3];
rz(2.7622278) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.9374342) q[2];
sx q[2];
rz(-1.5896475) q[2];
sx q[2];
rz(2.3834174) q[2];
rz(-0.6289064) q[3];
sx q[3];
rz(-0.40142504) q[3];
sx q[3];
rz(1.1531856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4085098) q[0];
sx q[0];
rz(-1.874431) q[0];
sx q[0];
rz(2.4531903) q[0];
rz(-0.06772659) q[1];
sx q[1];
rz(-1.7522782) q[1];
sx q[1];
rz(2.6115131) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2358658) q[0];
sx q[0];
rz(-1.4639963) q[0];
sx q[0];
rz(-1.279116) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7128503) q[2];
sx q[2];
rz(-1.0849761) q[2];
sx q[2];
rz(0.57810099) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.9024132) q[1];
sx q[1];
rz(-1.2352408) q[1];
sx q[1];
rz(-2.2526342) q[1];
x q[2];
rz(-1.9534555) q[3];
sx q[3];
rz(-0.2573765) q[3];
sx q[3];
rz(0.60636683) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.3337341) q[2];
sx q[2];
rz(-3.1299751) q[2];
sx q[2];
rz(-2.2401436) q[2];
rz(2.3060913) q[3];
sx q[3];
rz(-1.6136026) q[3];
sx q[3];
rz(-1.3114595) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19105844) q[0];
sx q[0];
rz(-1.5427417) q[0];
sx q[0];
rz(-0.78432551) q[0];
rz(-0.061231881) q[1];
sx q[1];
rz(-0.71413723) q[1];
sx q[1];
rz(3.004946) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2211321) q[0];
sx q[0];
rz(-1.4401299) q[0];
sx q[0];
rz(-1.0871825) q[0];
rz(2.6776671) q[2];
sx q[2];
rz(-1.8564965) q[2];
sx q[2];
rz(-2.6505016) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.5977051) q[1];
sx q[1];
rz(-1.502431) q[1];
sx q[1];
rz(-0.023936546) q[1];
x q[2];
rz(1.4284381) q[3];
sx q[3];
rz(-0.58066237) q[3];
sx q[3];
rz(-2.2992087) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1293929) q[2];
sx q[2];
rz(-2.2012074) q[2];
sx q[2];
rz(-0.56048918) q[2];
rz(0.012332049) q[3];
sx q[3];
rz(-0.90356946) q[3];
sx q[3];
rz(-1.0906609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
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
rz(2.7085768) q[0];
sx q[0];
rz(-2.5362159) q[0];
sx q[0];
rz(2.3204455) q[0];
rz(0.87617809) q[1];
sx q[1];
rz(-2.2416302) q[1];
sx q[1];
rz(1.4076153) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9579904) q[0];
sx q[0];
rz(-1.2215658) q[0];
sx q[0];
rz(2.7244085) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.9887092) q[2];
sx q[2];
rz(-1.313813) q[2];
sx q[2];
rz(-0.35494057) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.40286139) q[1];
sx q[1];
rz(-1.0051454) q[1];
sx q[1];
rz(2.9823142) q[1];
x q[2];
rz(-1.0293343) q[3];
sx q[3];
rz(-2.7280305) q[3];
sx q[3];
rz(-0.13973164) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.451482) q[2];
sx q[2];
rz(-1.2146981) q[2];
sx q[2];
rz(3.0991128) q[2];
rz(-2.5111607) q[3];
sx q[3];
rz(-0.63215956) q[3];
sx q[3];
rz(2.650034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7522488) q[0];
sx q[0];
rz(-1.2279953) q[0];
sx q[0];
rz(2.65843) q[0];
rz(1.0522316) q[1];
sx q[1];
rz(-1.9960884) q[1];
sx q[1];
rz(2.5767456) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4120867) q[0];
sx q[0];
rz(-0.54486638) q[0];
sx q[0];
rz(1.8338404) q[0];
x q[1];
rz(-3.0658423) q[2];
sx q[2];
rz(-1.9905914) q[2];
sx q[2];
rz(2.9043353) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.5833252) q[1];
sx q[1];
rz(-0.43154432) q[1];
sx q[1];
rz(1.7130909) q[1];
rz(-pi) q[2];
rz(-2.9103043) q[3];
sx q[3];
rz(-1.4196718) q[3];
sx q[3];
rz(-0.95201991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.57006449) q[2];
sx q[2];
rz(-1.0674942) q[2];
sx q[2];
rz(1.8035536) q[2];
rz(1.8367052) q[3];
sx q[3];
rz(-1.0675425) q[3];
sx q[3];
rz(-0.6750955) q[3];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9682482) q[0];
sx q[0];
rz(-1.4302379) q[0];
sx q[0];
rz(2.5937953) q[0];
rz(-2.3563747) q[1];
sx q[1];
rz(-1.3350057) q[1];
sx q[1];
rz(-0.26842591) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7086605) q[0];
sx q[0];
rz(-1.6527358) q[0];
sx q[0];
rz(2.8952778) q[0];
rz(-pi) q[1];
rz(2.3556049) q[2];
sx q[2];
rz(-0.95324264) q[2];
sx q[2];
rz(-2.3812889) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.0722326) q[1];
sx q[1];
rz(-1.6688804) q[1];
sx q[1];
rz(2.1588438) q[1];
rz(-2.8397452) q[3];
sx q[3];
rz(-2.0904623) q[3];
sx q[3];
rz(0.93320751) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.61775529) q[2];
sx q[2];
rz(-0.80344168) q[2];
sx q[2];
rz(-0.8141554) q[2];
rz(-2.7653149) q[3];
sx q[3];
rz(-1.9777931) q[3];
sx q[3];
rz(3.0686839) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58586621) q[0];
sx q[0];
rz(-1.7365475) q[0];
sx q[0];
rz(-2.1222173) q[0];
rz(0.85340071) q[1];
sx q[1];
rz(-1.1420206) q[1];
sx q[1];
rz(-0.44874915) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1837346) q[0];
sx q[0];
rz(-2.0439889) q[0];
sx q[0];
rz(0.26259043) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2737695) q[2];
sx q[2];
rz(-2.7542369) q[2];
sx q[2];
rz(-2.5572436) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.4510348) q[1];
sx q[1];
rz(-2.096855) q[1];
sx q[1];
rz(1.2760389) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3677424) q[3];
sx q[3];
rz(-1.0283096) q[3];
sx q[3];
rz(-2.0427996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.86924187) q[2];
sx q[2];
rz(-1.3808455) q[2];
sx q[2];
rz(0.31420079) q[2];
rz(-2.3172486) q[3];
sx q[3];
rz(-0.45212513) q[3];
sx q[3];
rz(-0.79469386) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.070351275) q[0];
sx q[0];
rz(-3.0817139) q[0];
sx q[0];
rz(-1.2605793) q[0];
rz(-0.64385995) q[1];
sx q[1];
rz(-1.2327797) q[1];
sx q[1];
rz(-0.018928122) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8277055) q[0];
sx q[0];
rz(-0.2642309) q[0];
sx q[0];
rz(3.0731191) q[0];
rz(2.5817714) q[2];
sx q[2];
rz(-0.18351843) q[2];
sx q[2];
rz(0.57932094) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.0156118) q[1];
sx q[1];
rz(-0.19128448) q[1];
sx q[1];
rz(0.26856883) q[1];
x q[2];
rz(-1.8065679) q[3];
sx q[3];
rz(-1.8528432) q[3];
sx q[3];
rz(2.3896133) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.5921322) q[2];
sx q[2];
rz(-0.38828725) q[2];
sx q[2];
rz(2.4712759) q[2];
rz(-0.51236764) q[3];
sx q[3];
rz(-1.3918326) q[3];
sx q[3];
rz(-1.7211154) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5877514) q[0];
sx q[0];
rz(-1.8974263) q[0];
sx q[0];
rz(0.49945369) q[0];
rz(-1.5669426) q[1];
sx q[1];
rz(-2.8630239) q[1];
sx q[1];
rz(1.0826899) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4944045) q[0];
sx q[0];
rz(-1.6196961) q[0];
sx q[0];
rz(-2.5542269) q[0];
rz(-2.322299) q[2];
sx q[2];
rz(-1.5466006) q[2];
sx q[2];
rz(-1.0011315) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-3.06304) q[1];
sx q[1];
rz(-1.5895956) q[1];
sx q[1];
rz(-2.5979554) q[1];
rz(-pi) q[2];
x q[2];
rz(0.70268481) q[3];
sx q[3];
rz(-2.6583238) q[3];
sx q[3];
rz(2.0600704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.7426976) q[2];
sx q[2];
rz(-1.164914) q[2];
sx q[2];
rz(-2.6297074) q[2];
rz(0.39294696) q[3];
sx q[3];
rz(-1.7397375) q[3];
sx q[3];
rz(-1.0569364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95505161) q[0];
sx q[0];
rz(-1.3165836) q[0];
sx q[0];
rz(-2.5008428) q[0];
rz(-2.4004249) q[1];
sx q[1];
rz(-0.82294958) q[1];
sx q[1];
rz(-0.23946147) q[1];
rz(0.95460931) q[2];
sx q[2];
rz(-1.2381427) q[2];
sx q[2];
rz(2.8340813) q[2];
rz(-0.74647222) q[3];
sx q[3];
rz(-1.6975879) q[3];
sx q[3];
rz(-3.0293037) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
