OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.0661434) q[0];
sx q[0];
rz(-2.0976522) q[0];
sx q[0];
rz(-0.010308417) q[0];
rz(-0.97996867) q[1];
sx q[1];
rz(4.5865321) q[1];
sx q[1];
rz(10.001339) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6714455) q[0];
sx q[0];
rz(-1.8785485) q[0];
sx q[0];
rz(-1.6977915) q[0];
rz(-pi) q[1];
rz(-0.89171191) q[2];
sx q[2];
rz(-2.0713965) q[2];
sx q[2];
rz(0.5989738) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.1548646) q[1];
sx q[1];
rz(-1.0467593) q[1];
sx q[1];
rz(-1.3681332) q[1];
rz(1.7553842) q[3];
sx q[3];
rz(-2.3999676) q[3];
sx q[3];
rz(-0.0024992873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.6171241) q[2];
sx q[2];
rz(-1.6813797) q[2];
sx q[2];
rz(1.2855533) q[2];
rz(-1.5517976) q[3];
sx q[3];
rz(-2.0701305) q[3];
sx q[3];
rz(-2.0901399) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0153506) q[0];
sx q[0];
rz(-1.9168357) q[0];
sx q[0];
rz(2.8958564) q[0];
rz(-2.083678) q[1];
sx q[1];
rz(-1.825288) q[1];
sx q[1];
rz(0.33448514) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4870105) q[0];
sx q[0];
rz(-1.0791057) q[0];
sx q[0];
rz(-1.0500748) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9742786) q[2];
sx q[2];
rz(-0.85131391) q[2];
sx q[2];
rz(3.1297504) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.3547438) q[1];
sx q[1];
rz(-1.9450099) q[1];
sx q[1];
rz(-1.1208865) q[1];
x q[2];
rz(1.2859551) q[3];
sx q[3];
rz(-2.1495616) q[3];
sx q[3];
rz(-2.9841686) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1841396) q[2];
sx q[2];
rz(-0.89749557) q[2];
sx q[2];
rz(-1.3857566) q[2];
rz(-0.98006788) q[3];
sx q[3];
rz(-2.6762784) q[3];
sx q[3];
rz(0.050962713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4362713) q[0];
sx q[0];
rz(-0.956981) q[0];
sx q[0];
rz(1.7514239) q[0];
rz(-0.35274371) q[1];
sx q[1];
rz(-2.2309062) q[1];
sx q[1];
rz(2.5211451) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1245956) q[0];
sx q[0];
rz(-0.11717883) q[0];
sx q[0];
rz(-2.4524053) q[0];
rz(-pi) q[1];
rz(-3.0376126) q[2];
sx q[2];
rz(-0.98538387) q[2];
sx q[2];
rz(2.3217161) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.31252334) q[1];
sx q[1];
rz(-1.8644973) q[1];
sx q[1];
rz(-1.2306661) q[1];
rz(-pi) q[2];
rz(0.87832344) q[3];
sx q[3];
rz(-1.9832423) q[3];
sx q[3];
rz(-0.45468047) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.1227526) q[2];
sx q[2];
rz(-2.7285125) q[2];
sx q[2];
rz(-1.8943465) q[2];
rz(2.9774169) q[3];
sx q[3];
rz(-1.434606) q[3];
sx q[3];
rz(2.5119787) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4267047) q[0];
sx q[0];
rz(-0.96788228) q[0];
sx q[0];
rz(-1.2870652) q[0];
rz(0.60028589) q[1];
sx q[1];
rz(-1.3515819) q[1];
sx q[1];
rz(-0.90369019) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0125192) q[0];
sx q[0];
rz(-0.12618574) q[0];
sx q[0];
rz(0.55380765) q[0];
rz(-pi) q[1];
rz(2.972347) q[2];
sx q[2];
rz(-2.7566559) q[2];
sx q[2];
rz(-2.9724309) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.43634448) q[1];
sx q[1];
rz(-1.7503993) q[1];
sx q[1];
rz(0.257538) q[1];
rz(-pi) q[2];
rz(-0.65151188) q[3];
sx q[3];
rz(-0.50262302) q[3];
sx q[3];
rz(0.16829695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.51957447) q[2];
sx q[2];
rz(-0.833424) q[2];
sx q[2];
rz(0.77345094) q[2];
rz(-1.2889688) q[3];
sx q[3];
rz(-0.52923146) q[3];
sx q[3];
rz(2.4066063) q[3];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2499823) q[0];
sx q[0];
rz(-2.1365428) q[0];
sx q[0];
rz(0.49215677) q[0];
rz(-0.53897578) q[1];
sx q[1];
rz(-1.69918) q[1];
sx q[1];
rz(2.0416562) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66456074) q[0];
sx q[0];
rz(-0.43433055) q[0];
sx q[0];
rz(-0.6566027) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.1308168) q[2];
sx q[2];
rz(-0.95913619) q[2];
sx q[2];
rz(-1.7623368) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.392994) q[1];
sx q[1];
rz(-1.6637319) q[1];
sx q[1];
rz(-0.8108906) q[1];
x q[2];
rz(1.2420916) q[3];
sx q[3];
rz(-2.0315758) q[3];
sx q[3];
rz(-2.8976909) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.8010572) q[2];
sx q[2];
rz(-0.33581442) q[2];
sx q[2];
rz(2.578793) q[2];
rz(-2.4417012) q[3];
sx q[3];
rz(-1.8507345) q[3];
sx q[3];
rz(-0.25115299) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3444779) q[0];
sx q[0];
rz(-2.4176702) q[0];
sx q[0];
rz(2.7864454) q[0];
rz(1.2190602) q[1];
sx q[1];
rz(-1.9025758) q[1];
sx q[1];
rz(-0.452279) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41846965) q[0];
sx q[0];
rz(-2.4112066) q[0];
sx q[0];
rz(0.95157051) q[0];
x q[1];
rz(-1.1638648) q[2];
sx q[2];
rz(-1.2135047) q[2];
sx q[2];
rz(2.1483764) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.86231316) q[1];
sx q[1];
rz(-1.8932027) q[1];
sx q[1];
rz(-3.1161948) q[1];
x q[2];
rz(2.5127453) q[3];
sx q[3];
rz(-1.4553242) q[3];
sx q[3];
rz(0.94146282) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.63048116) q[2];
sx q[2];
rz(-0.52383542) q[2];
sx q[2];
rz(3.0044978) q[2];
rz(-1.5824205) q[3];
sx q[3];
rz(-1.1538006) q[3];
sx q[3];
rz(-2.3649575) q[3];
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
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.95098507) q[0];
sx q[0];
rz(-2.3793716) q[0];
sx q[0];
rz(-1.5451587) q[0];
rz(-2.6243788) q[1];
sx q[1];
rz(-1.1511753) q[1];
sx q[1];
rz(0.98701611) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68970976) q[0];
sx q[0];
rz(-1.5919884) q[0];
sx q[0];
rz(-0.13084335) q[0];
x q[1];
rz(-2.8124078) q[2];
sx q[2];
rz(-0.98887695) q[2];
sx q[2];
rz(-2.0792368) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.92168857) q[1];
sx q[1];
rz(-2.7877586) q[1];
sx q[1];
rz(-2.4262731) q[1];
rz(-pi) q[2];
rz(1.889398) q[3];
sx q[3];
rz(-1.650839) q[3];
sx q[3];
rz(-2.6032676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.3551657) q[2];
sx q[2];
rz(-1.0309018) q[2];
sx q[2];
rz(-3.1332341) q[2];
rz(0.23056325) q[3];
sx q[3];
rz(-1.9356666) q[3];
sx q[3];
rz(1.6647388) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-0.01345988) q[0];
sx q[0];
rz(-0.020543329) q[0];
sx q[0];
rz(-1.531456) q[0];
rz(1.5229335) q[1];
sx q[1];
rz(-1.5956655) q[1];
sx q[1];
rz(1.5787554) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42199907) q[0];
sx q[0];
rz(-2.5223456) q[0];
sx q[0];
rz(0.25281711) q[0];
rz(0.35308102) q[2];
sx q[2];
rz(-0.39829474) q[2];
sx q[2];
rz(2.9495267) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.85962379) q[1];
sx q[1];
rz(-1.109237) q[1];
sx q[1];
rz(2.0412316) q[1];
rz(-pi) q[2];
rz(1.2501026) q[3];
sx q[3];
rz(-2.2134292) q[3];
sx q[3];
rz(-1.2144517) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.83850399) q[2];
sx q[2];
rz(-1.9946626) q[2];
sx q[2];
rz(-0.037503555) q[2];
rz(-2.904902) q[3];
sx q[3];
rz(-0.37875566) q[3];
sx q[3];
rz(-1.8481351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-2.8732052) q[0];
sx q[0];
rz(-0.78998843) q[0];
sx q[0];
rz(-1.0733806) q[0];
rz(2.7997596) q[1];
sx q[1];
rz(-0.57917246) q[1];
sx q[1];
rz(-2.5198708) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9699696) q[0];
sx q[0];
rz(-1.3116273) q[0];
sx q[0];
rz(-0.39879946) q[0];
rz(0.15839496) q[2];
sx q[2];
rz(-0.6662874) q[2];
sx q[2];
rz(1.0384384) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.3860064) q[1];
sx q[1];
rz(-1.2666128) q[1];
sx q[1];
rz(2.6261397) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.14071669) q[3];
sx q[3];
rz(-1.6575812) q[3];
sx q[3];
rz(-2.0279036) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.5198034) q[2];
sx q[2];
rz(-0.18814627) q[2];
sx q[2];
rz(-0.56358799) q[2];
rz(-1.311519) q[3];
sx q[3];
rz(-1.2389641) q[3];
sx q[3];
rz(0.81609503) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1821063) q[0];
sx q[0];
rz(-0.86903787) q[0];
sx q[0];
rz(2.2667789) q[0];
rz(-1.4080217) q[1];
sx q[1];
rz(-2.4751016) q[1];
sx q[1];
rz(-0.63546884) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68671526) q[0];
sx q[0];
rz(-1.2736763) q[0];
sx q[0];
rz(0.77154205) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9866139) q[2];
sx q[2];
rz(-1.8228616) q[2];
sx q[2];
rz(-2.3769555) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6136352) q[1];
sx q[1];
rz(-1.0571169) q[1];
sx q[1];
rz(0.46544816) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3395679) q[3];
sx q[3];
rz(-2.4200508) q[3];
sx q[3];
rz(-0.82593838) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.1982939) q[2];
sx q[2];
rz(-2.0512927) q[2];
sx q[2];
rz(2.3401006) q[2];
rz(1.9258026) q[3];
sx q[3];
rz(-0.91167584) q[3];
sx q[3];
rz(3.099814) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
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
rz(0.53032482) q[0];
sx q[0];
rz(-1.7014736) q[0];
sx q[0];
rz(0.3076719) q[0];
rz(1.8511741) q[1];
sx q[1];
rz(-2.1967874) q[1];
sx q[1];
rz(-2.8312942) q[1];
rz(0.18193131) q[2];
sx q[2];
rz(-1.7040737) q[2];
sx q[2];
rz(-2.2752442) q[2];
rz(-2.484816) q[3];
sx q[3];
rz(-1.4639502) q[3];
sx q[3];
rz(-1.7456036) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
