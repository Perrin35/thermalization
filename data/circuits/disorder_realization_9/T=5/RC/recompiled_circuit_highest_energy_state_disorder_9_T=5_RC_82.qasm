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
rz(3.4540662) q[1];
sx q[1];
rz(11.087505) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85017289) q[0];
sx q[0];
rz(-1.2787191) q[0];
sx q[0];
rz(-0.34698457) q[0];
x q[1];
rz(-1.0530627) q[2];
sx q[2];
rz(-1.3446756) q[2];
sx q[2];
rz(0.88498164) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.42740187) q[1];
sx q[1];
rz(-1.8334249) q[1];
sx q[1];
rz(-2.3841969) q[1];
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
rz(-2.5375157) q[2];
sx q[2];
rz(-0.61415577) q[2];
sx q[2];
rz(2.3094731) q[2];
rz(2.4845947) q[3];
sx q[3];
rz(-1.4703581) q[3];
sx q[3];
rz(2.7895797) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.326139) q[0];
sx q[0];
rz(-2.5322545) q[0];
sx q[0];
rz(-0.64999181) q[0];
rz(-3.0187712) q[1];
sx q[1];
rz(-0.42418066) q[1];
sx q[1];
rz(2.4688683) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8685659) q[0];
sx q[0];
rz(-1.582524) q[0];
sx q[0];
rz(-1.5290029) q[0];
x q[1];
rz(-0.71750516) q[2];
sx q[2];
rz(-0.40552545) q[2];
sx q[2];
rz(-2.6971291) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.0716463) q[1];
sx q[1];
rz(-2.5499747) q[1];
sx q[1];
rz(-1.2143192) q[1];
rz(-1.5244739) q[3];
sx q[3];
rz(-2.6213548) q[3];
sx q[3];
rz(-0.6689451) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.0279205) q[2];
sx q[2];
rz(-1.547926) q[2];
sx q[2];
rz(-0.38635722) q[2];
rz(-1.9733285) q[3];
sx q[3];
rz(-1.9100274) q[3];
sx q[3];
rz(-2.0333576) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89762178) q[0];
sx q[0];
rz(-1.6672927) q[0];
sx q[0];
rz(-0.058187159) q[0];
rz(2.8367786) q[1];
sx q[1];
rz(-2.0084281) q[1];
sx q[1];
rz(-2.286639) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8940272) q[0];
sx q[0];
rz(-1.6391048) q[0];
sx q[0];
rz(1.5724284) q[0];
x q[1];
rz(2.0570175) q[2];
sx q[2];
rz(-0.9409608) q[2];
sx q[2];
rz(-1.6597009) q[2];
x q[3];
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
x q[2];
rz(2.8775887) q[3];
sx q[3];
rz(-2.1101293) q[3];
sx q[3];
rz(-1.7738938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.052224) q[2];
sx q[2];
rz(-1.3244018) q[2];
sx q[2];
rz(-1.8734056) q[2];
rz(-2.7005699) q[3];
sx q[3];
rz(-2.5206168) q[3];
sx q[3];
rz(-2.2997901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12255254) q[0];
sx q[0];
rz(-0.95128107) q[0];
sx q[0];
rz(2.4265491) q[0];
rz(-2.7252281) q[1];
sx q[1];
rz(-2.6519897) q[1];
sx q[1];
rz(2.8474836) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.607632) q[0];
sx q[0];
rz(-2.2658969) q[0];
sx q[0];
rz(-3.0661723) q[0];
x q[1];
rz(1.3662228) q[2];
sx q[2];
rz(-1.6343717) q[2];
sx q[2];
rz(-0.83067375) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.4353883) q[1];
sx q[1];
rz(-1.7559373) q[1];
sx q[1];
rz(1.7035083) q[1];
rz(-0.58737602) q[3];
sx q[3];
rz(-2.0677635) q[3];
sx q[3];
rz(-2.9413829) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.9010345) q[2];
sx q[2];
rz(-1.3996539) q[2];
sx q[2];
rz(2.412839) q[2];
rz(0.48480222) q[3];
sx q[3];
rz(-1.0032283) q[3];
sx q[3];
rz(-0.17759594) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6840376) q[0];
sx q[0];
rz(-1.3223248) q[0];
sx q[0];
rz(0.27873248) q[0];
rz(3.0740671) q[1];
sx q[1];
rz(-2.4261256) q[1];
sx q[1];
rz(-2.6356437) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7161432) q[0];
sx q[0];
rz(-1.4977411) q[0];
sx q[0];
rz(0.55668932) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4691817) q[2];
sx q[2];
rz(-1.7676465) q[2];
sx q[2];
rz(0.54430279) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.5488779) q[1];
sx q[1];
rz(-1.028182) q[1];
sx q[1];
rz(2.1714669) q[1];
rz(-0.89754126) q[3];
sx q[3];
rz(-1.3607303) q[3];
sx q[3];
rz(2.1830934) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.2013596) q[2];
sx q[2];
rz(-1.4505922) q[2];
sx q[2];
rz(3.039321) q[2];
rz(0.30397948) q[3];
sx q[3];
rz(-2.961048) q[3];
sx q[3];
rz(2.6612018) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85222307) q[0];
sx q[0];
rz(-0.82813534) q[0];
sx q[0];
rz(2.4612259) q[0];
rz(-0.96087372) q[1];
sx q[1];
rz(-1.5545132) q[1];
sx q[1];
rz(-0.56112498) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8770804) q[0];
sx q[0];
rz(-1.3662369) q[0];
sx q[0];
rz(1.2007942) q[0];
rz(-pi) q[1];
x q[1];
rz(0.84662171) q[2];
sx q[2];
rz(-1.7270124) q[2];
sx q[2];
rz(-2.5515917) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.0099854) q[1];
sx q[1];
rz(-1.6254043) q[1];
sx q[1];
rz(-2.8427678) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4986491) q[3];
sx q[3];
rz(-0.74952945) q[3];
sx q[3];
rz(-1.3931605) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.23201021) q[2];
sx q[2];
rz(-2.1640615) q[2];
sx q[2];
rz(1.5781461) q[2];
rz(2.1252508) q[3];
sx q[3];
rz(-1.7137824) q[3];
sx q[3];
rz(2.0411172) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72961724) q[0];
sx q[0];
rz(-0.53891861) q[0];
sx q[0];
rz(2.6475661) q[0];
rz(-1.5771075) q[1];
sx q[1];
rz(-2.1421075) q[1];
sx q[1];
rz(-0.090242537) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0417953) q[0];
sx q[0];
rz(-1.3372412) q[0];
sx q[0];
rz(0.44333027) q[0];
x q[1];
rz(-0.97930564) q[2];
sx q[2];
rz(-2.3537209) q[2];
sx q[2];
rz(-0.56172127) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.75587326) q[1];
sx q[1];
rz(-1.8579646) q[1];
sx q[1];
rz(-1.1509291) q[1];
x q[2];
rz(-0.19321038) q[3];
sx q[3];
rz(-1.8487559) q[3];
sx q[3];
rz(1.4708468) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.6287441) q[2];
sx q[2];
rz(-2.5675842) q[2];
sx q[2];
rz(-2.4981456) q[2];
rz(-2.511034) q[3];
sx q[3];
rz(-0.75391155) q[3];
sx q[3];
rz(0.049230922) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.028320463) q[0];
sx q[0];
rz(-1.4406818) q[0];
sx q[0];
rz(1.2233618) q[0];
rz(-2.2471097) q[1];
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
rz(-1.2558738) q[0];
sx q[0];
rz(-0.71122716) q[0];
sx q[0];
rz(0.92899404) q[0];
x q[1];
rz(0.38739631) q[2];
sx q[2];
rz(-1.1715537) q[2];
sx q[2];
rz(0.42457595) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.8967197) q[1];
sx q[1];
rz(-1.280652) q[1];
sx q[1];
rz(-1.822132) q[1];
rz(0.98553879) q[3];
sx q[3];
rz(-1.7904591) q[3];
sx q[3];
rz(-2.1287763) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.2724096) q[2];
sx q[2];
rz(-0.86239186) q[2];
sx q[2];
rz(-1.4746846) q[2];
rz(0.21371755) q[3];
sx q[3];
rz(-1.4054047) q[3];
sx q[3];
rz(2.2739482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84233442) q[0];
sx q[0];
rz(-0.61115757) q[0];
sx q[0];
rz(0.70511955) q[0];
rz(0.57535386) q[1];
sx q[1];
rz(-1.4082785) q[1];
sx q[1];
rz(-2.5720678) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3045159) q[0];
sx q[0];
rz(-1.6479371) q[0];
sx q[0];
rz(1.6521554) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2766383) q[2];
sx q[2];
rz(-1.0964956) q[2];
sx q[2];
rz(-2.2612342) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.7895487) q[1];
sx q[1];
rz(-0.96992271) q[1];
sx q[1];
rz(-2.4070441) q[1];
rz(-pi) q[2];
rz(1.2678746) q[3];
sx q[3];
rz(-1.1008796) q[3];
sx q[3];
rz(2.5055714) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.9465955) q[2];
sx q[2];
rz(-2.1985998) q[2];
sx q[2];
rz(0.97879624) q[2];
rz(-0.43092522) q[3];
sx q[3];
rz(-2.6543591) q[3];
sx q[3];
rz(-2.0144958) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8579213) q[0];
sx q[0];
rz(-3.1223174) q[0];
sx q[0];
rz(1.4625782) q[0];
rz(1.0390394) q[1];
sx q[1];
rz(-1.3835399) q[1];
sx q[1];
rz(1.9336112) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6504575) q[0];
sx q[0];
rz(-1.8723328) q[0];
sx q[0];
rz(1.1745321) q[0];
rz(-pi) q[1];
rz(1.6287491) q[2];
sx q[2];
rz(-2.6160512) q[2];
sx q[2];
rz(-1.4144104) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.12700272) q[1];
sx q[1];
rz(-0.63204884) q[1];
sx q[1];
rz(1.3227425) q[1];
x q[2];
rz(1.1753756) q[3];
sx q[3];
rz(-1.5081769) q[3];
sx q[3];
rz(0.6400125) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.98564467) q[2];
sx q[2];
rz(-1.0040823) q[2];
sx q[2];
rz(-1.8211625) q[2];
rz(-2.7909347) q[3];
sx q[3];
rz(-0.81736332) q[3];
sx q[3];
rz(2.5613274) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.88631267) q[0];
sx q[0];
rz(-1.9552312) q[0];
sx q[0];
rz(2.2807688) q[0];
rz(1.4844004) q[1];
sx q[1];
rz(-1.3927554) q[1];
sx q[1];
rz(1.4295255) q[1];
rz(-1.4996573) q[2];
sx q[2];
rz(-1.9491458) q[2];
sx q[2];
rz(-0.083935621) q[2];
rz(0.30050011) q[3];
sx q[3];
rz(-0.90682744) q[3];
sx q[3];
rz(-1.711267) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
