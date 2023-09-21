OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.0619573) q[0];
sx q[0];
rz(-0.24793967) q[0];
sx q[0];
rz(2.3214582) q[0];
rz(0.040854383) q[1];
sx q[1];
rz(3.9305384) q[1];
sx q[1];
rz(9.5126704) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1577507) q[0];
sx q[0];
rz(-1.112126) q[0];
sx q[0];
rz(1.0755324) q[0];
rz(-pi) q[1];
rz(2.232993) q[2];
sx q[2];
rz(-2.5052921) q[2];
sx q[2];
rz(-1.5696021) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.80884113) q[1];
sx q[1];
rz(-1.0665227) q[1];
sx q[1];
rz(1.5732952) q[1];
x q[2];
rz(-1.5929234) q[3];
sx q[3];
rz(-2.8828354) q[3];
sx q[3];
rz(0.72507897) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.1408046) q[2];
sx q[2];
rz(-1.6892097) q[2];
sx q[2];
rz(0.56935707) q[2];
rz(1.6128444) q[3];
sx q[3];
rz(-2.5053535) q[3];
sx q[3];
rz(-1.7830085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59812087) q[0];
sx q[0];
rz(-0.16508979) q[0];
sx q[0];
rz(-0.55364451) q[0];
rz(-1.9042632) q[1];
sx q[1];
rz(-1.3668704) q[1];
sx q[1];
rz(1.9083317) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.935826) q[0];
sx q[0];
rz(-2.2201949) q[0];
sx q[0];
rz(2.4539024) q[0];
x q[1];
rz(2.9821175) q[2];
sx q[2];
rz(-2.4607686) q[2];
sx q[2];
rz(2.8603539) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.0765637) q[1];
sx q[1];
rz(-0.20670465) q[1];
sx q[1];
rz(1.4935342) q[1];
rz(-pi) q[2];
rz(-1.8096829) q[3];
sx q[3];
rz(-0.91730648) q[3];
sx q[3];
rz(1.7042004) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.0040940293) q[2];
sx q[2];
rz(-1.4592905) q[2];
sx q[2];
rz(-0.0022350524) q[2];
rz(-0.8301174) q[3];
sx q[3];
rz(-2.3486962) q[3];
sx q[3];
rz(1.8921651) q[3];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24901351) q[0];
sx q[0];
rz(-2.5254288) q[0];
sx q[0];
rz(2.2900443) q[0];
rz(0.7710723) q[1];
sx q[1];
rz(-1.4665736) q[1];
sx q[1];
rz(-3.070389) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.93807968) q[0];
sx q[0];
rz(-0.31103125) q[0];
sx q[0];
rz(2.2631893) q[0];
x q[1];
rz(-1.115828) q[2];
sx q[2];
rz(-1.3340545) q[2];
sx q[2];
rz(0.67982212) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.0752807) q[1];
sx q[1];
rz(-1.4367141) q[1];
sx q[1];
rz(-2.4877497) q[1];
rz(-pi) q[2];
x q[2];
rz(0.72910587) q[3];
sx q[3];
rz(-0.90568542) q[3];
sx q[3];
rz(2.2774746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.3123793) q[2];
sx q[2];
rz(-0.92386121) q[2];
sx q[2];
rz(2.6181347) q[2];
rz(0.33189014) q[3];
sx q[3];
rz(-3.0811946) q[3];
sx q[3];
rz(-2.6383242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4098542) q[0];
sx q[0];
rz(-0.24895746) q[0];
sx q[0];
rz(-0.82558924) q[0];
rz(-1.9748953) q[1];
sx q[1];
rz(-0.86667934) q[1];
sx q[1];
rz(1.4473787) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71890812) q[0];
sx q[0];
rz(-1.8005162) q[0];
sx q[0];
rz(-2.9740888) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7491313) q[2];
sx q[2];
rz(-1.5618088) q[2];
sx q[2];
rz(1.0967147) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.73562685) q[1];
sx q[1];
rz(-1.4117129) q[1];
sx q[1];
rz(-0.87079485) q[1];
rz(2.7778266) q[3];
sx q[3];
rz(-1.3437265) q[3];
sx q[3];
rz(-2.098339) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.46431413) q[2];
sx q[2];
rz(-1.364664) q[2];
sx q[2];
rz(0.14492598) q[2];
rz(2.1302917) q[3];
sx q[3];
rz(-2.5144808) q[3];
sx q[3];
rz(-3.1269126) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
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
rz(-0.99575627) q[0];
sx q[0];
rz(-3.0881112) q[0];
sx q[0];
rz(2.4575535) q[0];
rz(-1.9513291) q[1];
sx q[1];
rz(-0.66527706) q[1];
sx q[1];
rz(-1.9285944) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0048464674) q[0];
sx q[0];
rz(-0.46127013) q[0];
sx q[0];
rz(2.524551) q[0];
x q[1];
rz(2.5118675) q[2];
sx q[2];
rz(-1.0169528) q[2];
sx q[2];
rz(-0.50309203) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.2781093) q[1];
sx q[1];
rz(-2.8591595) q[1];
sx q[1];
rz(-2.7062098) q[1];
rz(-pi) q[2];
rz(2.1608859) q[3];
sx q[3];
rz(-0.31168918) q[3];
sx q[3];
rz(-2.2463472) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.6325536) q[2];
sx q[2];
rz(-1.2712487) q[2];
sx q[2];
rz(-0.63009134) q[2];
rz(0.57224327) q[3];
sx q[3];
rz(-2.4961491) q[3];
sx q[3];
rz(0.8852638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0260139) q[0];
sx q[0];
rz(-0.84340874) q[0];
sx q[0];
rz(-0.28636006) q[0];
rz(-2.5419366) q[1];
sx q[1];
rz(-1.1279794) q[1];
sx q[1];
rz(2.5568331) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2858546) q[0];
sx q[0];
rz(-1.2585088) q[0];
sx q[0];
rz(2.5166442) q[0];
rz(-pi) q[1];
rz(-1.5969959) q[2];
sx q[2];
rz(-0.95653406) q[2];
sx q[2];
rz(-2.7363077) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.11440052) q[1];
sx q[1];
rz(-2.235734) q[1];
sx q[1];
rz(-0.36892885) q[1];
rz(1.5201969) q[3];
sx q[3];
rz(-1.7369324) q[3];
sx q[3];
rz(1.43464) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.9770603) q[2];
sx q[2];
rz(-0.88447276) q[2];
sx q[2];
rz(-1.5131081) q[2];
rz(0.96308723) q[3];
sx q[3];
rz(-1.77308) q[3];
sx q[3];
rz(-2.517038) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.088783711) q[0];
sx q[0];
rz(-1.7099021) q[0];
sx q[0];
rz(-2.5928296) q[0];
rz(-3.0896297) q[1];
sx q[1];
rz(-0.49182645) q[1];
sx q[1];
rz(0.60639492) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0604027) q[0];
sx q[0];
rz(-2.5572544) q[0];
sx q[0];
rz(2.0085745) q[0];
rz(-pi) q[1];
x q[1];
rz(0.90784448) q[2];
sx q[2];
rz(-1.2693451) q[2];
sx q[2];
rz(0.19336685) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.7402621) q[1];
sx q[1];
rz(-1.4583424) q[1];
sx q[1];
rz(-0.80660352) q[1];
x q[2];
rz(2.7860836) q[3];
sx q[3];
rz(-1.1843369) q[3];
sx q[3];
rz(0.10458065) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.2018044) q[2];
sx q[2];
rz(-2.1494892) q[2];
sx q[2];
rz(-0.50841224) q[2];
rz(1.1172179) q[3];
sx q[3];
rz(-0.17337392) q[3];
sx q[3];
rz(-1.4846876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8023119) q[0];
sx q[0];
rz(-2.7957714) q[0];
sx q[0];
rz(2.3876277) q[0];
rz(2.1915961) q[1];
sx q[1];
rz(-1.6597304) q[1];
sx q[1];
rz(-2.9139013) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65093016) q[0];
sx q[0];
rz(-2.7224446) q[0];
sx q[0];
rz(1.7463669) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6954306) q[2];
sx q[2];
rz(-1.9036306) q[2];
sx q[2];
rz(-1.7762426) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.7714872) q[1];
sx q[1];
rz(-2.5574554) q[1];
sx q[1];
rz(2.7258956) q[1];
x q[2];
rz(-3.0759096) q[3];
sx q[3];
rz(-0.55995299) q[3];
sx q[3];
rz(-2.5064859) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.0630539) q[2];
sx q[2];
rz(-2.9324014) q[2];
sx q[2];
rz(0.83855808) q[2];
rz(-2.7770384) q[3];
sx q[3];
rz(-1.7611327) q[3];
sx q[3];
rz(2.3013039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7446328) q[0];
sx q[0];
rz(-1.7663706) q[0];
sx q[0];
rz(-0.80378419) q[0];
rz(-2.0869758) q[1];
sx q[1];
rz(-0.63916731) q[1];
sx q[1];
rz(-1.9931591) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9371985) q[0];
sx q[0];
rz(-0.075852595) q[0];
sx q[0];
rz(-1.8914468) q[0];
x q[1];
rz(0.9903332) q[2];
sx q[2];
rz(-0.60539421) q[2];
sx q[2];
rz(-0.47547542) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.941274) q[1];
sx q[1];
rz(-2.4793844) q[1];
sx q[1];
rz(-2.3108285) q[1];
rz(-pi) q[2];
rz(-1.6059444) q[3];
sx q[3];
rz(-1.9225549) q[3];
sx q[3];
rz(-2.7621321) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.049872963) q[2];
sx q[2];
rz(-2.1240978) q[2];
sx q[2];
rz(-2.3633374) q[2];
rz(-0.19566472) q[3];
sx q[3];
rz(-1.0396495) q[3];
sx q[3];
rz(-1.0218609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8675999) q[0];
sx q[0];
rz(-0.99853981) q[0];
sx q[0];
rz(0.25319779) q[0];
rz(-0.7397488) q[1];
sx q[1];
rz(-0.7363798) q[1];
sx q[1];
rz(-0.69828066) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.88903504) q[0];
sx q[0];
rz(-1.6425942) q[0];
sx q[0];
rz(-0.30307146) q[0];
x q[1];
rz(3.1168078) q[2];
sx q[2];
rz(-1.4966655) q[2];
sx q[2];
rz(-0.63400808) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2472898) q[1];
sx q[1];
rz(-1.3548046) q[1];
sx q[1];
rz(-2.0386019) q[1];
x q[2];
rz(1.0032759) q[3];
sx q[3];
rz(-1.365981) q[3];
sx q[3];
rz(1.2712461) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.11761052) q[2];
sx q[2];
rz(-0.88279804) q[2];
sx q[2];
rz(2.3013766) q[2];
rz(-2.501287) q[3];
sx q[3];
rz(-2.2655723) q[3];
sx q[3];
rz(1.1673814) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2587851) q[0];
sx q[0];
rz(-2.2866645) q[0];
sx q[0];
rz(-1.4186161) q[0];
rz(3.1124658) q[1];
sx q[1];
rz(-3.0283785) q[1];
sx q[1];
rz(1.821847) q[1];
rz(2.5393456) q[2];
sx q[2];
rz(-1.8613653) q[2];
sx q[2];
rz(-2.0801434) q[2];
rz(-0.15016951) q[3];
sx q[3];
rz(-0.88610813) q[3];
sx q[3];
rz(-0.025972493) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];