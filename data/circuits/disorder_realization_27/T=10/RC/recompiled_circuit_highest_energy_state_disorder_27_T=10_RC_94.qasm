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
rz(2.8994695) q[0];
sx q[0];
rz(-2.4162633) q[0];
sx q[0];
rz(-0.90670937) q[0];
rz(2.6746305) q[1];
sx q[1];
rz(-0.90277201) q[1];
sx q[1];
rz(-2.8977107) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7142993) q[0];
sx q[0];
rz(-2.1903945) q[0];
sx q[0];
rz(2.0626953) q[0];
x q[1];
rz(2.4087231) q[2];
sx q[2];
rz(-0.83031619) q[2];
sx q[2];
rz(0.24093369) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.4441057) q[1];
sx q[1];
rz(-0.19679697) q[1];
sx q[1];
rz(-2.7641731) q[1];
x q[2];
rz(-1.3451869) q[3];
sx q[3];
rz(-2.3757031) q[3];
sx q[3];
rz(-1.2768859) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.4563518) q[2];
sx q[2];
rz(-2.6754003) q[2];
sx q[2];
rz(-1.0249798) q[2];
rz(-0.13925615) q[3];
sx q[3];
rz(-1.9459629) q[3];
sx q[3];
rz(0.16417424) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76992947) q[0];
sx q[0];
rz(-2.0491845) q[0];
sx q[0];
rz(-1.2051693) q[0];
rz(-1.6734164) q[1];
sx q[1];
rz(-1.2638777) q[1];
sx q[1];
rz(-1.42043) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.842671) q[0];
sx q[0];
rz(-2.3793829) q[0];
sx q[0];
rz(-1.803714) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0864079) q[2];
sx q[2];
rz(-1.8351438) q[2];
sx q[2];
rz(2.4072539) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.5788889) q[1];
sx q[1];
rz(-1.3486224) q[1];
sx q[1];
rz(1.9166757) q[1];
rz(-0.74052625) q[3];
sx q[3];
rz(-1.9364001) q[3];
sx q[3];
rz(1.7056364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.18985192) q[2];
sx q[2];
rz(-0.18988374) q[2];
sx q[2];
rz(-2.3404549) q[2];
rz(0.43870157) q[3];
sx q[3];
rz(-1.2480241) q[3];
sx q[3];
rz(1.1130822) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37757117) q[0];
sx q[0];
rz(-0.29733297) q[0];
sx q[0];
rz(-1.1391033) q[0];
rz(0.26690075) q[1];
sx q[1];
rz(-1.2511988) q[1];
sx q[1];
rz(0.36531726) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8784155) q[0];
sx q[0];
rz(-0.52367822) q[0];
sx q[0];
rz(0.19864638) q[0];
rz(0.56198012) q[2];
sx q[2];
rz(-1.0022853) q[2];
sx q[2];
rz(-1.7719444) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.74148948) q[1];
sx q[1];
rz(-1.5129397) q[1];
sx q[1];
rz(-1.2025546) q[1];
rz(-pi) q[2];
rz(2.5015058) q[3];
sx q[3];
rz(-2.1800506) q[3];
sx q[3];
rz(-2.9971788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.2306564) q[2];
sx q[2];
rz(-2.4012884) q[2];
sx q[2];
rz(0.12348565) q[2];
rz(0.89928818) q[3];
sx q[3];
rz(-1.4495918) q[3];
sx q[3];
rz(1.2353157) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3673636) q[0];
sx q[0];
rz(-1.4643359) q[0];
sx q[0];
rz(2.2520219) q[0];
rz(2.3956237) q[1];
sx q[1];
rz(-1.4624701) q[1];
sx q[1];
rz(-1.3847345) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.810122) q[0];
sx q[0];
rz(-1.6056013) q[0];
sx q[0];
rz(-1.6286544) q[0];
rz(-pi) q[1];
rz(2.4998274) q[2];
sx q[2];
rz(-1.7671529) q[2];
sx q[2];
rz(-1.3349443) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.1390723) q[1];
sx q[1];
rz(-2.5031282) q[1];
sx q[1];
rz(2.8248766) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0880726) q[3];
sx q[3];
rz(-2.60703) q[3];
sx q[3];
rz(-0.10343119) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.7437462) q[2];
sx q[2];
rz(-0.83796871) q[2];
sx q[2];
rz(-1.0154593) q[2];
rz(-0.022653496) q[3];
sx q[3];
rz(-1.3709143) q[3];
sx q[3];
rz(-0.13815752) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37462336) q[0];
sx q[0];
rz(-1.2767108) q[0];
sx q[0];
rz(2.936777) q[0];
rz(-0.21518937) q[1];
sx q[1];
rz(-0.22061017) q[1];
sx q[1];
rz(0.87517103) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3615205) q[0];
sx q[0];
rz(-0.339739) q[0];
sx q[0];
rz(1.1113412) q[0];
x q[1];
rz(1.4007447) q[2];
sx q[2];
rz(-2.399246) q[2];
sx q[2];
rz(-0.60630732) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.4103247) q[1];
sx q[1];
rz(-1.350623) q[1];
sx q[1];
rz(0.30558821) q[1];
rz(2.2017893) q[3];
sx q[3];
rz(-1.6865589) q[3];
sx q[3];
rz(-2.1950658) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.508226) q[2];
sx q[2];
rz(-0.88819155) q[2];
sx q[2];
rz(1.7249031) q[2];
rz(-2.21777) q[3];
sx q[3];
rz(-2.1638736) q[3];
sx q[3];
rz(1.1317322) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1205587) q[0];
sx q[0];
rz(-2.4036305) q[0];
sx q[0];
rz(-0.84838867) q[0];
rz(1.5658762) q[1];
sx q[1];
rz(-1.8137685) q[1];
sx q[1];
rz(0.94589344) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83505982) q[0];
sx q[0];
rz(-0.9047821) q[0];
sx q[0];
rz(1.205797) q[0];
rz(-pi) q[1];
rz(-0.20393238) q[2];
sx q[2];
rz(-2.3205767) q[2];
sx q[2];
rz(-1.7199788) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.9528563) q[1];
sx q[1];
rz(-2.1486001) q[1];
sx q[1];
rz(1.0580479) q[1];
rz(-1.7686144) q[3];
sx q[3];
rz(-2.0536011) q[3];
sx q[3];
rz(1.5316666) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.1762323) q[2];
sx q[2];
rz(-1.9051899) q[2];
sx q[2];
rz(1.4005093) q[2];
rz(1.4817574) q[3];
sx q[3];
rz(-1.7912495) q[3];
sx q[3];
rz(0.19937936) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7401212) q[0];
sx q[0];
rz(-0.49982247) q[0];
sx q[0];
rz(1.2116145) q[0];
rz(0.47677332) q[1];
sx q[1];
rz(-1.2766653) q[1];
sx q[1];
rz(-1.5062987) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20902987) q[0];
sx q[0];
rz(-1.2878006) q[0];
sx q[0];
rz(-2.1542633) q[0];
rz(2.2231323) q[2];
sx q[2];
rz(-2.7428584) q[2];
sx q[2];
rz(-0.14130558) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.16852779) q[1];
sx q[1];
rz(-2.2650411) q[1];
sx q[1];
rz(1.4100875) q[1];
x q[2];
rz(2.6166418) q[3];
sx q[3];
rz(-0.54407507) q[3];
sx q[3];
rz(-1.3348188) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3393042) q[2];
sx q[2];
rz(-1.5172639) q[2];
sx q[2];
rz(-2.9586672) q[2];
rz(0.66169345) q[3];
sx q[3];
rz(-1.2122093) q[3];
sx q[3];
rz(-0.70484149) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6905007) q[0];
sx q[0];
rz(-1.6766312) q[0];
sx q[0];
rz(-0.14725421) q[0];
rz(0.14689771) q[1];
sx q[1];
rz(-2.2203827) q[1];
sx q[1];
rz(-1.9619092) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3578276) q[0];
sx q[0];
rz(-1.234878) q[0];
sx q[0];
rz(-0.442105) q[0];
rz(-2.7034892) q[2];
sx q[2];
rz(-1.1313651) q[2];
sx q[2];
rz(2.8218486) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.6322286) q[1];
sx q[1];
rz(-1.6809789) q[1];
sx q[1];
rz(2.5446507) q[1];
rz(-pi) q[2];
rz(3.0249765) q[3];
sx q[3];
rz(-0.75748235) q[3];
sx q[3];
rz(2.1471786) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.94072056) q[2];
sx q[2];
rz(-2.2042037) q[2];
sx q[2];
rz(-0.76784039) q[2];
rz(0.52669865) q[3];
sx q[3];
rz(-2.0898762) q[3];
sx q[3];
rz(0.55829486) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.679477) q[0];
sx q[0];
rz(-0.18590346) q[0];
sx q[0];
rz(-2.0515077) q[0];
rz(-1.7915626) q[1];
sx q[1];
rz(-1.7410024) q[1];
sx q[1];
rz(2.7166691) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37575484) q[0];
sx q[0];
rz(-2.0415487) q[0];
sx q[0];
rz(-0.063675675) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0504136) q[2];
sx q[2];
rz(-2.489739) q[2];
sx q[2];
rz(-0.21287795) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.7560727) q[1];
sx q[1];
rz(-1.2282097) q[1];
sx q[1];
rz(-0.48348654) q[1];
rz(-0.70885371) q[3];
sx q[3];
rz(-1.9836622) q[3];
sx q[3];
rz(-0.23504892) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.2316042) q[2];
sx q[2];
rz(-0.51077545) q[2];
sx q[2];
rz(-2.6166022) q[2];
rz(-1.3548343) q[3];
sx q[3];
rz(-2.2460008) q[3];
sx q[3];
rz(1.3625328) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43596426) q[0];
sx q[0];
rz(-2.4810677) q[0];
sx q[0];
rz(-3.0349162) q[0];
rz(-1.0542997) q[1];
sx q[1];
rz(-1.432447) q[1];
sx q[1];
rz(-1.4073102) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76864734) q[0];
sx q[0];
rz(-2.1539186) q[0];
sx q[0];
rz(-2.5685422) q[0];
rz(-pi) q[1];
rz(2.5564747) q[2];
sx q[2];
rz(-0.74293929) q[2];
sx q[2];
rz(-1.9415426) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.52508721) q[1];
sx q[1];
rz(-0.27180755) q[1];
sx q[1];
rz(0.96801438) q[1];
rz(2.1945215) q[3];
sx q[3];
rz(-2.2184988) q[3];
sx q[3];
rz(0.050681678) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.6537031) q[2];
sx q[2];
rz(-2.6370878) q[2];
sx q[2];
rz(-2.8686236) q[2];
rz(-0.87131396) q[3];
sx q[3];
rz(-1.6815691) q[3];
sx q[3];
rz(-0.13450809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7731666) q[0];
sx q[0];
rz(-1.1933403) q[0];
sx q[0];
rz(-2.2884952) q[0];
rz(2.7586965) q[1];
sx q[1];
rz(-1.8538414) q[1];
sx q[1];
rz(-0.014898653) q[1];
rz(0.56131526) q[2];
sx q[2];
rz(-1.9519162) q[2];
sx q[2];
rz(-1.219195) q[2];
rz(3.0440108) q[3];
sx q[3];
rz(-2.3392936) q[3];
sx q[3];
rz(0.38264075) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
