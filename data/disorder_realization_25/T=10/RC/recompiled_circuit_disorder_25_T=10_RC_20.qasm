OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.0796354) q[0];
sx q[0];
rz(-2.893653) q[0];
sx q[0];
rz(-2.3214582) q[0];
rz(-3.1007383) q[1];
sx q[1];
rz(-0.78894579) q[1];
sx q[1];
rz(-0.087892428) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9633006) q[0];
sx q[0];
rz(-1.1305729) q[0];
sx q[0];
rz(-2.630169) q[0];
rz(2.232993) q[2];
sx q[2];
rz(-0.63630051) q[2];
sx q[2];
rz(-1.5719906) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.3808448) q[1];
sx q[1];
rz(-1.5686085) q[1];
sx q[1];
rz(-0.50427498) q[1];
rz(-pi) q[2];
rz(0.0058562972) q[3];
sx q[3];
rz(-1.3121038) q[3];
sx q[3];
rz(-0.70219016) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1408046) q[2];
sx q[2];
rz(-1.452383) q[2];
sx q[2];
rz(0.56935707) q[2];
rz(1.5287483) q[3];
sx q[3];
rz(-0.6362392) q[3];
sx q[3];
rz(1.3585842) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59812087) q[0];
sx q[0];
rz(-0.16508979) q[0];
sx q[0];
rz(-2.5879481) q[0];
rz(-1.2373295) q[1];
sx q[1];
rz(-1.3668704) q[1];
sx q[1];
rz(1.233261) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.935826) q[0];
sx q[0];
rz(-2.2201949) q[0];
sx q[0];
rz(-0.68769023) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.1594752) q[2];
sx q[2];
rz(-2.4607686) q[2];
sx q[2];
rz(2.8603539) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.0765637) q[1];
sx q[1];
rz(-0.20670465) q[1];
sx q[1];
rz(-1.4935342) q[1];
x q[2];
rz(0.66744653) q[3];
sx q[3];
rz(-1.3818041) q[3];
sx q[3];
rz(-2.8611956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(3.1374986) q[2];
sx q[2];
rz(-1.4592905) q[2];
sx q[2];
rz(0.0022350524) q[2];
rz(2.3114752) q[3];
sx q[3];
rz(-0.79289645) q[3];
sx q[3];
rz(1.2494276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8925791) q[0];
sx q[0];
rz(-0.61616388) q[0];
sx q[0];
rz(0.85154831) q[0];
rz(-0.7710723) q[1];
sx q[1];
rz(-1.4665736) q[1];
sx q[1];
rz(3.070389) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.203513) q[0];
sx q[0];
rz(-0.31103125) q[0];
sx q[0];
rz(2.2631893) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.26239563) q[2];
sx q[2];
rz(-2.0121644) q[2];
sx q[2];
rz(1.0052094) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.0752807) q[1];
sx q[1];
rz(-1.7048786) q[1];
sx q[1];
rz(2.4877497) q[1];
x q[2];
rz(0.72910587) q[3];
sx q[3];
rz(-2.2359072) q[3];
sx q[3];
rz(-2.2774746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.8292134) q[2];
sx q[2];
rz(-2.2177314) q[2];
sx q[2];
rz(0.52345792) q[2];
rz(2.8097025) q[3];
sx q[3];
rz(-3.0811946) q[3];
sx q[3];
rz(2.6383242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(1.7317384) q[0];
sx q[0];
rz(-2.8926352) q[0];
sx q[0];
rz(0.82558924) q[0];
rz(1.9748953) q[1];
sx q[1];
rz(-0.86667934) q[1];
sx q[1];
rz(-1.4473787) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3576413) q[0];
sx q[0];
rz(-0.28342993) q[0];
sx q[0];
rz(2.1901603) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6214192) q[2];
sx q[2];
rz(-2.9630337) q[2];
sx q[2];
rz(0.52390097) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.64926681) q[1];
sx q[1];
rz(-0.71486231) q[1];
sx q[1];
rz(1.8148755) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7778266) q[3];
sx q[3];
rz(-1.7978661) q[3];
sx q[3];
rz(-2.098339) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.6772785) q[2];
sx q[2];
rz(-1.7769287) q[2];
sx q[2];
rz(-0.14492598) q[2];
rz(-1.011301) q[3];
sx q[3];
rz(-2.5144808) q[3];
sx q[3];
rz(0.01468006) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1458364) q[0];
sx q[0];
rz(-0.053481426) q[0];
sx q[0];
rz(0.68403912) q[0];
rz(-1.9513291) q[1];
sx q[1];
rz(-0.66527706) q[1];
sx q[1];
rz(-1.9285944) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0048464674) q[0];
sx q[0];
rz(-0.46127013) q[0];
sx q[0];
rz(0.61704163) q[0];
rz(-pi) q[1];
x q[1];
rz(0.91765399) q[2];
sx q[2];
rz(-2.0954164) q[2];
sx q[2];
rz(1.4337002) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8634833) q[1];
sx q[1];
rz(-0.28243318) q[1];
sx q[1];
rz(-0.43538283) q[1];
rz(-pi) q[2];
rz(-2.1608859) q[3];
sx q[3];
rz(-2.8299035) q[3];
sx q[3];
rz(0.89524549) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.6325536) q[2];
sx q[2];
rz(-1.8703439) q[2];
sx q[2];
rz(0.63009134) q[2];
rz(0.57224327) q[3];
sx q[3];
rz(-0.64544353) q[3];
sx q[3];
rz(2.2563289) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
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
rz(1.1155788) q[0];
sx q[0];
rz(-0.84340874) q[0];
sx q[0];
rz(-2.8552326) q[0];
rz(-2.5419366) q[1];
sx q[1];
rz(-1.1279794) q[1];
sx q[1];
rz(-0.58475959) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8557381) q[0];
sx q[0];
rz(-1.2585088) q[0];
sx q[0];
rz(0.62494846) q[0];
x q[1];
rz(-1.5969959) q[2];
sx q[2];
rz(-0.95653406) q[2];
sx q[2];
rz(0.40528497) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.6905745) q[1];
sx q[1];
rz(-1.283053) q[1];
sx q[1];
rz(-0.87177999) q[1];
rz(-pi) q[2];
rz(-1.6213958) q[3];
sx q[3];
rz(-1.7369324) q[3];
sx q[3];
rz(-1.7069526) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9770603) q[2];
sx q[2];
rz(-0.88447276) q[2];
sx q[2];
rz(-1.5131081) q[2];
rz(0.96308723) q[3];
sx q[3];
rz(-1.3685127) q[3];
sx q[3];
rz(-0.62455463) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.088783711) q[0];
sx q[0];
rz(-1.4316906) q[0];
sx q[0];
rz(-0.54876304) q[0];
rz(0.051963003) q[1];
sx q[1];
rz(-2.6497662) q[1];
sx q[1];
rz(-0.60639492) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.279778) q[0];
sx q[0];
rz(-1.3347515) q[0];
sx q[0];
rz(1.0311014) q[0];
rz(-pi) q[1];
rz(2.2337482) q[2];
sx q[2];
rz(-1.8722476) q[2];
sx q[2];
rz(0.19336685) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.4013306) q[1];
sx q[1];
rz(-1.4583424) q[1];
sx q[1];
rz(2.3349891) q[1];
x q[2];
rz(-1.9803195) q[3];
sx q[3];
rz(-1.2425353) q[3];
sx q[3];
rz(1.3271774) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.2018044) q[2];
sx q[2];
rz(-0.99210343) q[2];
sx q[2];
rz(-0.50841224) q[2];
rz(2.0243747) q[3];
sx q[3];
rz(-0.17337392) q[3];
sx q[3];
rz(-1.6569051) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8023119) q[0];
sx q[0];
rz(-2.7957714) q[0];
sx q[0];
rz(0.75396496) q[0];
rz(0.94999653) q[1];
sx q[1];
rz(-1.6597304) q[1];
sx q[1];
rz(2.9139013) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.75922155) q[0];
sx q[0];
rz(-1.4996487) q[0];
sx q[0];
rz(1.1573777) q[0];
x q[1];
rz(2.4661602) q[2];
sx q[2];
rz(-2.5917412) q[2];
sx q[2];
rz(-2.7477802) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.5880809) q[1];
sx q[1];
rz(-1.3462102) q[1];
sx q[1];
rz(-0.54393804) q[1];
rz(0.55898198) q[3];
sx q[3];
rz(-1.6056656) q[3];
sx q[3];
rz(-2.1502286) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(3.0630539) q[2];
sx q[2];
rz(-0.20919122) q[2];
sx q[2];
rz(-0.83855808) q[2];
rz(2.7770384) q[3];
sx q[3];
rz(-1.7611327) q[3];
sx q[3];
rz(-2.3013039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39695981) q[0];
sx q[0];
rz(-1.3752221) q[0];
sx q[0];
rz(-0.80378419) q[0];
rz(1.0546168) q[1];
sx q[1];
rz(-2.5024253) q[1];
sx q[1];
rz(1.9931591) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20439416) q[0];
sx q[0];
rz(-3.0657401) q[0];
sx q[0];
rz(1.2501459) q[0];
x q[1];
rz(-0.36275136) q[2];
sx q[2];
rz(-2.0667549) q[2];
sx q[2];
rz(0.19778684) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.8879537) q[1];
sx q[1];
rz(-1.1432853) q[1];
sx q[1];
rz(1.0484139) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5356482) q[3];
sx q[3];
rz(-1.9225549) q[3];
sx q[3];
rz(2.7621321) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.049872963) q[2];
sx q[2];
rz(-1.0174948) q[2];
sx q[2];
rz(2.3633374) q[2];
rz(0.19566472) q[3];
sx q[3];
rz(-2.1019432) q[3];
sx q[3];
rz(2.1197317) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2739928) q[0];
sx q[0];
rz(-0.99853981) q[0];
sx q[0];
rz(-2.8883949) q[0];
rz(-2.4018438) q[1];
sx q[1];
rz(-2.4052129) q[1];
sx q[1];
rz(2.443312) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90726501) q[0];
sx q[0];
rz(-0.31120473) q[0];
sx q[0];
rz(-2.9051203) q[0];
rz(-pi) q[1];
rz(-1.8928705) q[2];
sx q[2];
rz(-0.078157166) q[2];
sx q[2];
rz(-2.8305778) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.56863927) q[1];
sx q[1];
rz(-1.1146953) q[1];
sx q[1];
rz(-2.9005463) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9000557) q[3];
sx q[3];
rz(-1.0165443) q[3];
sx q[3];
rz(0.17061558) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.11761052) q[2];
sx q[2];
rz(-0.88279804) q[2];
sx q[2];
rz(-0.84021604) q[2];
rz(2.501287) q[3];
sx q[3];
rz(-0.87602031) q[3];
sx q[3];
rz(1.1673814) q[3];
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
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8828076) q[0];
sx q[0];
rz(-2.2866645) q[0];
sx q[0];
rz(-1.4186161) q[0];
rz(-3.1124658) q[1];
sx q[1];
rz(-0.11321414) q[1];
sx q[1];
rz(-1.3197457) q[1];
rz(2.5393456) q[2];
sx q[2];
rz(-1.8613653) q[2];
sx q[2];
rz(-2.0801434) q[2];
rz(-2.9914231) q[3];
sx q[3];
rz(-2.2554845) q[3];
sx q[3];
rz(3.1156202) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
