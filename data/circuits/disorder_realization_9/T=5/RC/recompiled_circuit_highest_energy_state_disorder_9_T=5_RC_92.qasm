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
rz(0.63737386) q[0];
sx q[0];
rz(-2.0847991) q[0];
sx q[0];
rz(-0.771653) q[0];
rz(-0.15089384) q[1];
sx q[1];
rz(-0.3124736) q[1];
sx q[1];
rz(-1.6627275) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2914198) q[0];
sx q[0];
rz(-1.2787191) q[0];
sx q[0];
rz(-0.34698457) q[0];
rz(-pi) q[1];
x q[1];
rz(0.25881473) q[2];
sx q[2];
rz(-2.0740899) q[2];
sx q[2];
rz(-0.81282114) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.42740187) q[1];
sx q[1];
rz(-1.8334249) q[1];
sx q[1];
rz(2.3841969) q[1];
rz(-pi) q[2];
rz(1.7467878) q[3];
sx q[3];
rz(-0.8991836) q[3];
sx q[3];
rz(-1.0351537) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.60407698) q[2];
sx q[2];
rz(-2.5274369) q[2];
sx q[2];
rz(0.83211952) q[2];
rz(-0.65699792) q[3];
sx q[3];
rz(-1.6712345) q[3];
sx q[3];
rz(-2.7895797) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81545365) q[0];
sx q[0];
rz(-2.5322545) q[0];
sx q[0];
rz(-2.4916008) q[0];
rz(-0.12282148) q[1];
sx q[1];
rz(-0.42418066) q[1];
sx q[1];
rz(-2.4688683) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.024351748) q[0];
sx q[0];
rz(-3.0981859) q[0];
sx q[0];
rz(-1.2971334) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.31284903) q[2];
sx q[2];
rz(-1.8331844) q[2];
sx q[2];
rz(-2.69115) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.9422226) q[1];
sx q[1];
rz(-1.7666715) q[1];
sx q[1];
rz(2.1327399) q[1];
x q[2];
rz(-2.0905714) q[3];
sx q[3];
rz(-1.5938164) q[3];
sx q[3];
rz(2.1995403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.1136721) q[2];
sx q[2];
rz(-1.5936667) q[2];
sx q[2];
rz(-0.38635722) q[2];
rz(-1.1682642) q[3];
sx q[3];
rz(-1.2315653) q[3];
sx q[3];
rz(-2.0333576) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89762178) q[0];
sx q[0];
rz(-1.6672927) q[0];
sx q[0];
rz(3.0834055) q[0];
rz(0.30481401) q[1];
sx q[1];
rz(-2.0084281) q[1];
sx q[1];
rz(-0.85495368) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2475654) q[0];
sx q[0];
rz(-1.5024878) q[0];
sx q[0];
rz(-1.5691643) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0845752) q[2];
sx q[2];
rz(-2.2006319) q[2];
sx q[2];
rz(1.6597009) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3173609) q[1];
sx q[1];
rz(-0.79272017) q[1];
sx q[1];
rz(0.2893682) q[1];
rz(-pi) q[2];
rz(1.0157973) q[3];
sx q[3];
rz(-1.3449752) q[3];
sx q[3];
rz(2.8005536) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.052224) q[2];
sx q[2];
rz(-1.3244018) q[2];
sx q[2];
rz(-1.8734056) q[2];
rz(-2.7005699) q[3];
sx q[3];
rz(-0.62097582) q[3];
sx q[3];
rz(2.2997901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.12255254) q[0];
sx q[0];
rz(-2.1903116) q[0];
sx q[0];
rz(-0.71504354) q[0];
rz(0.41636458) q[1];
sx q[1];
rz(-2.6519897) q[1];
sx q[1];
rz(2.8474836) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65140066) q[0];
sx q[0];
rz(-2.4430902) q[0];
sx q[0];
rz(1.480689) q[0];
x q[1];
rz(-1.2671183) q[2];
sx q[2];
rz(-0.21409245) q[2];
sx q[2];
rz(-2.10432) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.3015705) q[1];
sx q[1];
rz(-1.4403655) q[1];
sx q[1];
rz(-2.954847) q[1];
x q[2];
rz(0.58737602) q[3];
sx q[3];
rz(-2.0677635) q[3];
sx q[3];
rz(2.9413829) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.9010345) q[2];
sx q[2];
rz(-1.3996539) q[2];
sx q[2];
rz(-0.72875363) q[2];
rz(-0.48480222) q[3];
sx q[3];
rz(-1.0032283) q[3];
sx q[3];
rz(-2.9639967) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4575551) q[0];
sx q[0];
rz(-1.8192679) q[0];
sx q[0];
rz(2.8628602) q[0];
rz(0.067525603) q[1];
sx q[1];
rz(-0.7154671) q[1];
sx q[1];
rz(-2.6356437) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.099952817) q[0];
sx q[0];
rz(-1.0157662) q[0];
sx q[0];
rz(-1.4848079) q[0];
rz(2.9437482) q[2];
sx q[2];
rz(-1.4711507) q[2];
sx q[2];
rz(2.0951592) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.7796353) q[1];
sx q[1];
rz(-2.07603) q[1];
sx q[1];
rz(-0.63118668) q[1];
rz(0.89754126) q[3];
sx q[3];
rz(-1.3607303) q[3];
sx q[3];
rz(-2.1830934) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.94023306) q[2];
sx q[2];
rz(-1.6910005) q[2];
sx q[2];
rz(-0.10227164) q[2];
rz(-0.30397948) q[3];
sx q[3];
rz(-0.18054466) q[3];
sx q[3];
rz(2.6612018) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2893696) q[0];
sx q[0];
rz(-0.82813534) q[0];
sx q[0];
rz(-0.68036675) q[0];
rz(-0.96087372) q[1];
sx q[1];
rz(-1.5545132) q[1];
sx q[1];
rz(-0.56112498) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38491098) q[0];
sx q[0];
rz(-1.2088684) q[0];
sx q[0];
rz(0.2189526) q[0];
rz(-1.3374001) q[2];
sx q[2];
rz(-0.73783685) q[2];
sx q[2];
rz(-0.8066752) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.7192156) q[1];
sx q[1];
rz(-1.869162) q[1];
sx q[1];
rz(1.6279312) q[1];
rz(-pi) q[2];
rz(-2.5013728) q[3];
sx q[3];
rz(-1.9915765) q[3];
sx q[3];
rz(-2.8176378) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.23201021) q[2];
sx q[2];
rz(-0.97753111) q[2];
sx q[2];
rz(-1.5781461) q[2];
rz(1.0163418) q[3];
sx q[3];
rz(-1.4278102) q[3];
sx q[3];
rz(2.0411172) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72961724) q[0];
sx q[0];
rz(-0.53891861) q[0];
sx q[0];
rz(0.49402657) q[0];
rz(-1.5644851) q[1];
sx q[1];
rz(-0.99948519) q[1];
sx q[1];
rz(-0.090242537) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2170885) q[0];
sx q[0];
rz(-2.6441433) q[0];
sx q[0];
rz(0.50636943) q[0];
rz(-pi) q[1];
rz(2.162287) q[2];
sx q[2];
rz(-0.78787178) q[2];
sx q[2];
rz(0.56172127) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.7612865) q[1];
sx q[1];
rz(-0.50380987) q[1];
sx q[1];
rz(0.94382091) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9483823) q[3];
sx q[3];
rz(-1.2928367) q[3];
sx q[3];
rz(-1.4708468) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.51284853) q[2];
sx q[2];
rz(-2.5675842) q[2];
sx q[2];
rz(-2.4981456) q[2];
rz(-0.63055864) q[3];
sx q[3];
rz(-2.3876811) q[3];
sx q[3];
rz(-3.0923617) q[3];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1132722) q[0];
sx q[0];
rz(-1.4406818) q[0];
sx q[0];
rz(-1.2233618) q[0];
rz(0.89448294) q[1];
sx q[1];
rz(-0.62791413) q[1];
sx q[1];
rz(1.1462513) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.664331) q[0];
sx q[0];
rz(-1.0205752) q[0];
sx q[0];
rz(2.6653637) q[0];
x q[1];
rz(2.3010761) q[2];
sx q[2];
rz(-0.54893755) q[2];
sx q[2];
rz(-2.7565589) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.97630298) q[1];
sx q[1];
rz(-0.38152757) q[1];
sx q[1];
rz(2.4470727) q[1];
rz(-pi) q[2];
rz(-1.9548863) q[3];
sx q[3];
rz(-2.5210104) q[3];
sx q[3];
rz(2.9012802) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.2724096) q[2];
sx q[2];
rz(-0.86239186) q[2];
sx q[2];
rz(1.4746846) q[2];
rz(0.21371755) q[3];
sx q[3];
rz(-1.4054047) q[3];
sx q[3];
rz(2.2739482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84233442) q[0];
sx q[0];
rz(-0.61115757) q[0];
sx q[0];
rz(0.70511955) q[0];
rz(-2.5662388) q[1];
sx q[1];
rz(-1.4082785) q[1];
sx q[1];
rz(-2.5720678) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0234719) q[0];
sx q[0];
rz(-3.0295353) q[0];
sx q[0];
rz(-0.81046354) q[0];
rz(2.6492278) q[2];
sx q[2];
rz(-1.8316744) q[2];
sx q[2];
rz(-2.3136669) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.77806738) q[1];
sx q[1];
rz(-0.91202867) q[1];
sx q[1];
rz(-2.3450065) q[1];
x q[2];
rz(1.2678746) q[3];
sx q[3];
rz(-1.1008796) q[3];
sx q[3];
rz(2.5055714) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.19499714) q[2];
sx q[2];
rz(-0.94299287) q[2];
sx q[2];
rz(-0.97879624) q[2];
rz(2.7106674) q[3];
sx q[3];
rz(-0.48723358) q[3];
sx q[3];
rz(2.0144958) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8579213) q[0];
sx q[0];
rz(-0.019275276) q[0];
sx q[0];
rz(-1.6790144) q[0];
rz(1.0390394) q[1];
sx q[1];
rz(-1.7580527) q[1];
sx q[1];
rz(-1.9336112) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20328612) q[0];
sx q[0];
rz(-1.948256) q[0];
sx q[0];
rz(-0.32517821) q[0];
x q[1];
rz(-0.033577888) q[2];
sx q[2];
rz(-2.0953645) q[2];
sx q[2];
rz(-1.3474423) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.12700272) q[1];
sx q[1];
rz(-2.5095438) q[1];
sx q[1];
rz(1.8188502) q[1];
x q[2];
rz(0.067840067) q[3];
sx q[3];
rz(-1.1761936) q[3];
sx q[3];
rz(2.2369235) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.155948) q[2];
sx q[2];
rz(-2.1375103) q[2];
sx q[2];
rz(1.8211625) q[2];
rz(-2.7909347) q[3];
sx q[3];
rz(-0.81736332) q[3];
sx q[3];
rz(-0.58026522) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.25528) q[0];
sx q[0];
rz(-1.1863615) q[0];
sx q[0];
rz(-0.86082389) q[0];
rz(1.6571922) q[1];
sx q[1];
rz(-1.7488372) q[1];
sx q[1];
rz(-1.7120672) q[1];
rz(-0.17694494) q[2];
sx q[2];
rz(-2.7569303) q[2];
sx q[2];
rz(-3.0349572) q[2];
rz(-0.30050011) q[3];
sx q[3];
rz(-2.2347652) q[3];
sx q[3];
rz(1.4303257) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
