OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.3818504) q[0];
sx q[0];
rz(2.3072825) q[0];
sx q[0];
rz(6.351525) q[0];
rz(-0.66032687) q[1];
sx q[1];
rz(-0.84815174) q[1];
sx q[1];
rz(0.037820427) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67536394) q[0];
sx q[0];
rz(-0.20875202) q[0];
sx q[0];
rz(-1.1842313) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.90934609) q[2];
sx q[2];
rz(-1.6644018) q[2];
sx q[2];
rz(-2.0441165) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.0953857) q[1];
sx q[1];
rz(-2.3933105) q[1];
sx q[1];
rz(-3.066643) q[1];
rz(-1.2245523) q[3];
sx q[3];
rz(-1.4284705) q[3];
sx q[3];
rz(3.0179253) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.35090703) q[2];
sx q[2];
rz(-2.6280792) q[2];
sx q[2];
rz(-1.2834056) q[2];
rz(3.0154199) q[3];
sx q[3];
rz(-1.4161371) q[3];
sx q[3];
rz(-0.090601966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.693817) q[0];
sx q[0];
rz(-0.87647804) q[0];
sx q[0];
rz(-0.59666657) q[0];
rz(1.5860575) q[1];
sx q[1];
rz(-1.3417473) q[1];
sx q[1];
rz(-1.3635051) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1024433) q[0];
sx q[0];
rz(-1.1457232) q[0];
sx q[0];
rz(2.5545679) q[0];
rz(-pi) q[1];
rz(-1.6215789) q[2];
sx q[2];
rz(-2.1261458) q[2];
sx q[2];
rz(0.83088779) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.1494257) q[1];
sx q[1];
rz(-0.70632315) q[1];
sx q[1];
rz(-0.43744932) q[1];
rz(-pi) q[2];
rz(0.84888208) q[3];
sx q[3];
rz(-0.97182453) q[3];
sx q[3];
rz(0.90585432) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.8186701) q[2];
sx q[2];
rz(-2.1408036) q[2];
sx q[2];
rz(-2.4070516) q[2];
rz(-2.5143886) q[3];
sx q[3];
rz(-1.6278798) q[3];
sx q[3];
rz(-2.396615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4194141) q[0];
sx q[0];
rz(-0.83291554) q[0];
sx q[0];
rz(-0.27221361) q[0];
rz(2.294337) q[1];
sx q[1];
rz(-1.3191185) q[1];
sx q[1];
rz(0.31262696) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74041849) q[0];
sx q[0];
rz(-2.9156988) q[0];
sx q[0];
rz(-0.24143879) q[0];
x q[1];
rz(-0.69487822) q[2];
sx q[2];
rz(-2.4701397) q[2];
sx q[2];
rz(1.7365255) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.4362267) q[1];
sx q[1];
rz(-1.6402906) q[1];
sx q[1];
rz(1.8314929) q[1];
rz(-0.39194312) q[3];
sx q[3];
rz(-0.80229811) q[3];
sx q[3];
rz(-2.9585569) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.12038885) q[2];
sx q[2];
rz(-0.52336064) q[2];
sx q[2];
rz(2.3266501) q[2];
rz(-0.96873823) q[3];
sx q[3];
rz(-1.0453984) q[3];
sx q[3];
rz(-0.43280861) q[3];
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
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75354904) q[0];
sx q[0];
rz(-0.20064813) q[0];
sx q[0];
rz(-0.62336212) q[0];
rz(2.3240044) q[1];
sx q[1];
rz(-2.8249884) q[1];
sx q[1];
rz(-0.049291074) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3885956) q[0];
sx q[0];
rz(-1.6191282) q[0];
sx q[0];
rz(2.3301793) q[0];
rz(2.8231499) q[2];
sx q[2];
rz(-0.88469425) q[2];
sx q[2];
rz(2.7903914) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4781487) q[1];
sx q[1];
rz(-1.2423007) q[1];
sx q[1];
rz(-2.8847787) q[1];
x q[2];
rz(0.87818273) q[3];
sx q[3];
rz(-2.4947824) q[3];
sx q[3];
rz(-1.1991024) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.7307044) q[2];
sx q[2];
rz(-1.5521908) q[2];
sx q[2];
rz(0.98497406) q[2];
rz(-0.90304053) q[3];
sx q[3];
rz(-1.1629546) q[3];
sx q[3];
rz(-0.27339098) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(-0.34489283) q[0];
sx q[0];
rz(-1.6051689) q[0];
sx q[0];
rz(3.0539736) q[0];
rz(2.9852729) q[1];
sx q[1];
rz(-0.55115288) q[1];
sx q[1];
rz(0.87096754) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1317469) q[0];
sx q[0];
rz(-1.5994659) q[0];
sx q[0];
rz(-2.8542551) q[0];
rz(-1.9268553) q[2];
sx q[2];
rz(-2.7241754) q[2];
sx q[2];
rz(-2.2603214) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.37967967) q[1];
sx q[1];
rz(-1.2907791) q[1];
sx q[1];
rz(0.67358394) q[1];
x q[2];
rz(1.3887651) q[3];
sx q[3];
rz(-2.3689299) q[3];
sx q[3];
rz(2.3327737) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.0481723) q[2];
sx q[2];
rz(-0.8478567) q[2];
sx q[2];
rz(-0.35287738) q[2];
rz(2.2757163) q[3];
sx q[3];
rz(-2.4398949) q[3];
sx q[3];
rz(-2.849259) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46309328) q[0];
sx q[0];
rz(-2.0149639) q[0];
sx q[0];
rz(0.10678664) q[0];
rz(-1.1865901) q[1];
sx q[1];
rz(-2.1148966) q[1];
sx q[1];
rz(-2.1170763) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9232641) q[0];
sx q[0];
rz(-2.5946293) q[0];
sx q[0];
rz(-0.57759072) q[0];
rz(-1.1149939) q[2];
sx q[2];
rz(-2.1264646) q[2];
sx q[2];
rz(1.9089886) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.7926327) q[1];
sx q[1];
rz(-1.0981202) q[1];
sx q[1];
rz(2.2036168) q[1];
x q[2];
rz(2.5173353) q[3];
sx q[3];
rz(-0.58371021) q[3];
sx q[3];
rz(0.39238413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.391905) q[2];
sx q[2];
rz(-1.2140032) q[2];
sx q[2];
rz(-0.8708896) q[2];
rz(-1.8093367) q[3];
sx q[3];
rz(-0.94227666) q[3];
sx q[3];
rz(-1.4107305) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1720599) q[0];
sx q[0];
rz(-1.4071858) q[0];
sx q[0];
rz(-2.5906738) q[0];
rz(-2.6761966) q[1];
sx q[1];
rz(-0.67133633) q[1];
sx q[1];
rz(-0.23682061) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.119334) q[0];
sx q[0];
rz(-0.10611457) q[0];
sx q[0];
rz(-1.2073713) q[0];
rz(-1.4332438) q[2];
sx q[2];
rz(-1.7195065) q[2];
sx q[2];
rz(-2.6580722) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.9501818) q[1];
sx q[1];
rz(-1.4958053) q[1];
sx q[1];
rz(-1.300315) q[1];
x q[2];
rz(-1.9067494) q[3];
sx q[3];
rz(-0.7630322) q[3];
sx q[3];
rz(1.2515765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6992496) q[2];
sx q[2];
rz(-2.4599059) q[2];
sx q[2];
rz(-2.701475) q[2];
rz(-1.0951428) q[3];
sx q[3];
rz(-1.4705855) q[3];
sx q[3];
rz(1.346689) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8286164) q[0];
sx q[0];
rz(-1.799311) q[0];
sx q[0];
rz(3.112088) q[0];
rz(0.94738952) q[1];
sx q[1];
rz(-0.82059971) q[1];
sx q[1];
rz(-0.11880076) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4122075) q[0];
sx q[0];
rz(-2.8280624) q[0];
sx q[0];
rz(3.0010812) q[0];
rz(-pi) q[1];
rz(2.115909) q[2];
sx q[2];
rz(-1.9695749) q[2];
sx q[2];
rz(1.1023956) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5549705) q[1];
sx q[1];
rz(-1.8645617) q[1];
sx q[1];
rz(-0.48420669) q[1];
rz(-pi) q[2];
rz(1.0766255) q[3];
sx q[3];
rz(-1.4255187) q[3];
sx q[3];
rz(0.52500099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.7416731) q[2];
sx q[2];
rz(-2.6779149) q[2];
sx q[2];
rz(1.5734394) q[2];
rz(-1.167477) q[3];
sx q[3];
rz(-1.7220595) q[3];
sx q[3];
rz(-1.8673816) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90299273) q[0];
sx q[0];
rz(-2.3165343) q[0];
sx q[0];
rz(2.9883244) q[0];
rz(2.080147) q[1];
sx q[1];
rz(-2.5669284) q[1];
sx q[1];
rz(-1.9877888) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.783219) q[0];
sx q[0];
rz(-1.8405387) q[0];
sx q[0];
rz(-2.7620035) q[0];
rz(-pi) q[1];
rz(-2.597528) q[2];
sx q[2];
rz(-0.96368921) q[2];
sx q[2];
rz(2.0008759) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.8968618) q[1];
sx q[1];
rz(-1.3618999) q[1];
sx q[1];
rz(3.0848068) q[1];
x q[2];
rz(2.2581402) q[3];
sx q[3];
rz(-1.6646107) q[3];
sx q[3];
rz(2.1283619) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.29685059) q[2];
sx q[2];
rz(-1.1721609) q[2];
sx q[2];
rz(1.5273757) q[2];
rz(-0.99689233) q[3];
sx q[3];
rz(-1.6485873) q[3];
sx q[3];
rz(-2.2629288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24755724) q[0];
sx q[0];
rz(-1.0972801) q[0];
sx q[0];
rz(-2.9515008) q[0];
rz(-2.4709573) q[1];
sx q[1];
rz(-1.9452483) q[1];
sx q[1];
rz(1.6533096) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7877486) q[0];
sx q[0];
rz(-1.6962546) q[0];
sx q[0];
rz(-1.5397443) q[0];
rz(-0.44234862) q[2];
sx q[2];
rz(-1.8162236) q[2];
sx q[2];
rz(-0.52877141) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.59939811) q[1];
sx q[1];
rz(-2.290526) q[1];
sx q[1];
rz(-0.69516121) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2404352) q[3];
sx q[3];
rz(-0.81448758) q[3];
sx q[3];
rz(0.27064532) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.0037447475) q[2];
sx q[2];
rz(-2.2403084) q[2];
sx q[2];
rz(0.89912644) q[2];
rz(-2.6265465) q[3];
sx q[3];
rz(-0.47596541) q[3];
sx q[3];
rz(2.9639444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0777733) q[0];
sx q[0];
rz(-1.780092) q[0];
sx q[0];
rz(-0.5583981) q[0];
rz(2.8521815) q[1];
sx q[1];
rz(-2.2402973) q[1];
sx q[1];
rz(-1.4351861) q[1];
rz(0.32439705) q[2];
sx q[2];
rz(-0.73642052) q[2];
sx q[2];
rz(-3.0036075) q[2];
rz(-1.1424941) q[3];
sx q[3];
rz(-2.8430568) q[3];
sx q[3];
rz(-1.6381016) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
