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
rz(-2.5042188) q[0];
sx q[0];
rz(-1.0567935) q[0];
sx q[0];
rz(-2.3699397) q[0];
rz(2.9906988) q[1];
sx q[1];
rz(-2.8291191) q[1];
sx q[1];
rz(1.6627275) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7485086) q[0];
sx q[0];
rz(-2.6918852) q[0];
sx q[0];
rz(0.72400399) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0530627) q[2];
sx q[2];
rz(-1.7969171) q[2];
sx q[2];
rz(0.88498164) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.7141908) q[1];
sx q[1];
rz(-1.3081677) q[1];
sx q[1];
rz(2.3841969) q[1];
x q[2];
rz(-2.9247901) q[3];
sx q[3];
rz(-2.4507782) q[3];
sx q[3];
rz(1.8280713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.5375157) q[2];
sx q[2];
rz(-2.5274369) q[2];
sx q[2];
rz(2.3094731) q[2];
rz(-0.65699792) q[3];
sx q[3];
rz(-1.4703581) q[3];
sx q[3];
rz(-0.35201296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.326139) q[0];
sx q[0];
rz(-0.6093381) q[0];
sx q[0];
rz(-2.4916008) q[0];
rz(3.0187712) q[1];
sx q[1];
rz(-2.717412) q[1];
sx q[1];
rz(2.4688683) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2730267) q[0];
sx q[0];
rz(-1.5590686) q[0];
sx q[0];
rz(1.6125897) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8459199) q[2];
sx q[2];
rz(-1.2689991) q[2];
sx q[2];
rz(-1.2040621) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.4933659) q[1];
sx q[1];
rz(-1.0208482) q[1];
sx q[1];
rz(-0.2303161) q[1];
rz(2.0905714) q[3];
sx q[3];
rz(-1.5477763) q[3];
sx q[3];
rz(2.1995403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.0279205) q[2];
sx q[2];
rz(-1.547926) q[2];
sx q[2];
rz(-0.38635722) q[2];
rz(1.9733285) q[3];
sx q[3];
rz(-1.9100274) q[3];
sx q[3];
rz(2.0333576) q[3];
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
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89762178) q[0];
sx q[0];
rz(-1.4742999) q[0];
sx q[0];
rz(3.0834055) q[0];
rz(-2.8367786) q[1];
sx q[1];
rz(-1.1331646) q[1];
sx q[1];
rz(-2.286639) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2475654) q[0];
sx q[0];
rz(-1.5024878) q[0];
sx q[0];
rz(-1.5724284) q[0];
rz(-pi) q[1];
rz(-2.0570175) q[2];
sx q[2];
rz(-0.9409608) q[2];
sx q[2];
rz(-1.4818918) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.1890167) q[1];
sx q[1];
rz(-1.3661279) q[1];
sx q[1];
rz(0.7714899) q[1];
x q[2];
rz(1.1596572) q[3];
sx q[3];
rz(-2.5469031) q[3];
sx q[3];
rz(-0.88319544) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.0893687) q[2];
sx q[2];
rz(-1.3244018) q[2];
sx q[2];
rz(-1.8734056) q[2];
rz(2.7005699) q[3];
sx q[3];
rz(-0.62097582) q[3];
sx q[3];
rz(0.84180251) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0190401) q[0];
sx q[0];
rz(-2.1903116) q[0];
sx q[0];
rz(2.4265491) q[0];
rz(2.7252281) q[1];
sx q[1];
rz(-0.48960296) q[1];
sx q[1];
rz(-0.29410902) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.607632) q[0];
sx q[0];
rz(-0.87569571) q[0];
sx q[0];
rz(3.0661723) q[0];
rz(-pi) q[1];
x q[1];
rz(0.064925504) q[2];
sx q[2];
rz(-1.366642) q[2];
sx q[2];
rz(2.4146508) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.84002217) q[1];
sx q[1];
rz(-1.4403655) q[1];
sx q[1];
rz(2.954847) q[1];
x q[2];
rz(0.77463051) q[3];
sx q[3];
rz(-0.74997682) q[3];
sx q[3];
rz(2.3924654) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.2405582) q[2];
sx q[2];
rz(-1.7419387) q[2];
sx q[2];
rz(-2.412839) q[2];
rz(-2.6567904) q[3];
sx q[3];
rz(-2.1383643) q[3];
sx q[3];
rz(0.17759594) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4575551) q[0];
sx q[0];
rz(-1.8192679) q[0];
sx q[0];
rz(0.27873248) q[0];
rz(3.0740671) q[1];
sx q[1];
rz(-2.4261256) q[1];
sx q[1];
rz(0.50594893) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4254494) q[0];
sx q[0];
rz(-1.4977411) q[0];
sx q[0];
rz(0.55668932) q[0];
rz(-0.19784446) q[2];
sx q[2];
rz(-1.6704419) q[2];
sx q[2];
rz(-2.0951592) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.3619574) q[1];
sx q[1];
rz(-2.07603) q[1];
sx q[1];
rz(-2.510406) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9002731) q[3];
sx q[3];
rz(-2.4412324) q[3];
sx q[3];
rz(-0.86802378) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.2013596) q[2];
sx q[2];
rz(-1.6910005) q[2];
sx q[2];
rz(-0.10227164) q[2];
rz(-0.30397948) q[3];
sx q[3];
rz(-2.961048) q[3];
sx q[3];
rz(-2.6612018) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85222307) q[0];
sx q[0];
rz(-0.82813534) q[0];
sx q[0];
rz(0.68036675) q[0];
rz(0.96087372) q[1];
sx q[1];
rz(-1.5870794) q[1];
sx q[1];
rz(-0.56112498) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8770804) q[0];
sx q[0];
rz(-1.7753557) q[0];
sx q[0];
rz(-1.9407985) q[0];
x q[1];
rz(2.9343453) q[2];
sx q[2];
rz(-2.2842513) q[2];
sx q[2];
rz(-2.0240473) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.0099854) q[1];
sx q[1];
rz(-1.5161884) q[1];
sx q[1];
rz(-2.8427678) q[1];
rz(-pi) q[2];
rz(-2.0797727) q[3];
sx q[3];
rz(-0.99409311) q[3];
sx q[3];
rz(0.95147607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9095824) q[2];
sx q[2];
rz(-2.1640615) q[2];
sx q[2];
rz(1.5781461) q[2];
rz(2.1252508) q[3];
sx q[3];
rz(-1.4278102) q[3];
sx q[3];
rz(1.1004755) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-0.72961724) q[0];
sx q[0];
rz(-2.602674) q[0];
sx q[0];
rz(2.6475661) q[0];
rz(1.5771075) q[1];
sx q[1];
rz(-2.1421075) q[1];
sx q[1];
rz(0.090242537) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3615389) q[0];
sx q[0];
rz(-1.1403196) q[0];
sx q[0];
rz(-1.3132877) q[0];
rz(-pi) q[1];
rz(0.51076575) q[2];
sx q[2];
rz(-0.94168451) q[2];
sx q[2];
rz(0.19925403) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.7612865) q[1];
sx q[1];
rz(-2.6377828) q[1];
sx q[1];
rz(-0.94382091) q[1];
x q[2];
rz(-2.163104) q[3];
sx q[3];
rz(-2.8045125) q[3];
sx q[3];
rz(-0.85142985) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.51284853) q[2];
sx q[2];
rz(-0.57400846) q[2];
sx q[2];
rz(-2.4981456) q[2];
rz(2.511034) q[3];
sx q[3];
rz(-0.75391155) q[3];
sx q[3];
rz(3.0923617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.028320463) q[0];
sx q[0];
rz(-1.7009108) q[0];
sx q[0];
rz(1.9182308) q[0];
rz(2.2471097) q[1];
sx q[1];
rz(-0.62791413) q[1];
sx q[1];
rz(1.9953413) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.664331) q[0];
sx q[0];
rz(-2.1210175) q[0];
sx q[0];
rz(-0.47622891) q[0];
rz(2.3010761) q[2];
sx q[2];
rz(-2.5926551) q[2];
sx q[2];
rz(-0.38503371) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.24487296) q[1];
sx q[1];
rz(-1.280652) q[1];
sx q[1];
rz(-1.822132) q[1];
x q[2];
rz(-2.1560539) q[3];
sx q[3];
rz(-1.7904591) q[3];
sx q[3];
rz(-2.1287763) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.869183) q[2];
sx q[2];
rz(-2.2792008) q[2];
sx q[2];
rz(-1.4746846) q[2];
rz(-2.9278751) q[3];
sx q[3];
rz(-1.4054047) q[3];
sx q[3];
rz(-0.86764446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84233442) q[0];
sx q[0];
rz(-0.61115757) q[0];
sx q[0];
rz(-0.70511955) q[0];
rz(-2.5662388) q[1];
sx q[1];
rz(-1.4082785) q[1];
sx q[1];
rz(-2.5720678) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8370767) q[0];
sx q[0];
rz(-1.4936556) q[0];
sx q[0];
rz(-1.4894372) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6274983) q[2];
sx q[2];
rz(-2.5894458) q[2];
sx q[2];
rz(-0.29468003) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.25334206) q[1];
sx q[1];
rz(-0.98501316) q[1];
sx q[1];
rz(-0.82510494) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2678746) q[3];
sx q[3];
rz(-2.0407131) q[3];
sx q[3];
rz(-0.63602122) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.9465955) q[2];
sx q[2];
rz(-0.94299287) q[2];
sx q[2];
rz(-0.97879624) q[2];
rz(2.7106674) q[3];
sx q[3];
rz(-2.6543591) q[3];
sx q[3];
rz(-2.0144958) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8579213) q[0];
sx q[0];
rz(-3.1223174) q[0];
sx q[0];
rz(1.6790144) q[0];
rz(-2.1025533) q[1];
sx q[1];
rz(-1.7580527) q[1];
sx q[1];
rz(-1.9336112) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6504575) q[0];
sx q[0];
rz(-1.8723328) q[0];
sx q[0];
rz(-1.9670606) q[0];
x q[1];
rz(3.1080148) q[2];
sx q[2];
rz(-1.0462282) q[2];
sx q[2];
rz(1.3474423) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.9644248) q[1];
sx q[1];
rz(-0.96098122) q[1];
sx q[1];
rz(-0.17788203) q[1];
x q[2];
rz(-1.7321587) q[3];
sx q[3];
rz(-0.40009016) q[3];
sx q[3];
rz(2.0619947) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.98564467) q[2];
sx q[2];
rz(-2.1375103) q[2];
sx q[2];
rz(1.3204302) q[2];
rz(-2.7909347) q[3];
sx q[3];
rz(-2.3242293) q[3];
sx q[3];
rz(-2.5613274) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(-0.88631267) q[0];
sx q[0];
rz(-1.1863615) q[0];
sx q[0];
rz(-0.86082389) q[0];
rz(1.6571922) q[1];
sx q[1];
rz(-1.7488372) q[1];
sx q[1];
rz(-1.7120672) q[1];
rz(-0.37921956) q[2];
sx q[2];
rz(-1.6368964) q[2];
sx q[2];
rz(1.5131769) q[2];
rz(2.8410925) q[3];
sx q[3];
rz(-2.2347652) q[3];
sx q[3];
rz(1.4303257) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
