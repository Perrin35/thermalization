OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.7597423) q[0];
sx q[0];
rz(-2.3072825) q[0];
sx q[0];
rz(0.068339737) q[0];
rz(2.4812658) q[1];
sx q[1];
rz(3.9897444) q[1];
sx q[1];
rz(6.2453649) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8605182) q[0];
sx q[0];
rz(-1.3776508) q[0];
sx q[0];
rz(-3.0618969) q[0];
rz(-2.2322466) q[2];
sx q[2];
rz(-1.4771909) q[2];
sx q[2];
rz(1.0974761) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.046207) q[1];
sx q[1];
rz(-0.74828212) q[1];
sx q[1];
rz(3.066643) q[1];
rz(-1.171265) q[3];
sx q[3];
rz(-0.37326187) q[3];
sx q[3];
rz(-1.319862) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.35090703) q[2];
sx q[2];
rz(-2.6280792) q[2];
sx q[2];
rz(1.8581871) q[2];
rz(3.0154199) q[3];
sx q[3];
rz(-1.4161371) q[3];
sx q[3];
rz(3.0509907) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44777563) q[0];
sx q[0];
rz(-0.87647804) q[0];
sx q[0];
rz(-2.5449261) q[0];
rz(1.5860575) q[1];
sx q[1];
rz(-1.7998453) q[1];
sx q[1];
rz(-1.7780875) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1024433) q[0];
sx q[0];
rz(-1.1457232) q[0];
sx q[0];
rz(-0.58702472) q[0];
x q[1];
rz(-1.6215789) q[2];
sx q[2];
rz(-1.0154468) q[2];
sx q[2];
rz(-0.83088779) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9047782) q[1];
sx q[1];
rz(-1.292255) q[1];
sx q[1];
rz(-2.4836471) q[1];
rz(-pi) q[2];
rz(-2.2927106) q[3];
sx q[3];
rz(-2.1697681) q[3];
sx q[3];
rz(2.2357383) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.8186701) q[2];
sx q[2];
rz(-2.1408036) q[2];
sx q[2];
rz(-2.4070516) q[2];
rz(-0.62720403) q[3];
sx q[3];
rz(-1.5137129) q[3];
sx q[3];
rz(0.74497765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4194141) q[0];
sx q[0];
rz(-0.83291554) q[0];
sx q[0];
rz(0.27221361) q[0];
rz(0.84725562) q[1];
sx q[1];
rz(-1.8224742) q[1];
sx q[1];
rz(0.31262696) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74041849) q[0];
sx q[0];
rz(-2.9156988) q[0];
sx q[0];
rz(0.24143879) q[0];
rz(-pi) q[1];
rz(0.54801268) q[2];
sx q[2];
rz(-1.9805038) q[2];
sx q[2];
rz(0.41248413) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.12007512) q[1];
sx q[1];
rz(-2.8719963) q[1];
sx q[1];
rz(-1.3070379) q[1];
rz(1.1945046) q[3];
sx q[3];
rz(-0.84405758) q[3];
sx q[3];
rz(2.4220667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(3.0212038) q[2];
sx q[2];
rz(-2.618232) q[2];
sx q[2];
rz(0.81494251) q[2];
rz(-2.1728544) q[3];
sx q[3];
rz(-1.0453984) q[3];
sx q[3];
rz(0.43280861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3880436) q[0];
sx q[0];
rz(-2.9409445) q[0];
sx q[0];
rz(-0.62336212) q[0];
rz(0.81758824) q[1];
sx q[1];
rz(-0.31660429) q[1];
sx q[1];
rz(3.0923016) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.374644) q[0];
sx q[0];
rz(-2.3809803) q[0];
sx q[0];
rz(1.6409372) q[0];
x q[1];
rz(2.8231499) q[2];
sx q[2];
rz(-2.2568984) q[2];
sx q[2];
rz(0.35120121) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(3.1218159) q[1];
sx q[1];
rz(-0.41408086) q[1];
sx q[1];
rz(-2.2112234) q[1];
x q[2];
rz(2.0972581) q[3];
sx q[3];
rz(-1.1757848) q[3];
sx q[3];
rz(2.928424) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.7307044) q[2];
sx q[2];
rz(-1.5894019) q[2];
sx q[2];
rz(-2.1566186) q[2];
rz(-0.90304053) q[3];
sx q[3];
rz(-1.9786381) q[3];
sx q[3];
rz(0.27339098) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7966998) q[0];
sx q[0];
rz(-1.6051689) q[0];
sx q[0];
rz(-0.087619089) q[0];
rz(-2.9852729) q[1];
sx q[1];
rz(-0.55115288) q[1];
sx q[1];
rz(2.2706251) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.464251) q[0];
sx q[0];
rz(-2.8528677) q[0];
sx q[0];
rz(3.0407453) q[0];
rz(-pi) q[1];
rz(-1.9268553) q[2];
sx q[2];
rz(-0.41741727) q[2];
sx q[2];
rz(2.2603214) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.7333784) q[1];
sx q[1];
rz(-2.2136662) q[1];
sx q[1];
rz(1.9233568) q[1];
x q[2];
rz(-0.17467588) q[3];
sx q[3];
rz(-2.3274765) q[3];
sx q[3];
rz(-2.0811618) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.093420371) q[2];
sx q[2];
rz(-2.293736) q[2];
sx q[2];
rz(0.35287738) q[2];
rz(0.86587632) q[3];
sx q[3];
rz(-2.4398949) q[3];
sx q[3];
rz(2.849259) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
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
rz(-0.46309328) q[0];
sx q[0];
rz(-1.1266288) q[0];
sx q[0];
rz(-3.034806) q[0];
rz(-1.9550025) q[1];
sx q[1];
rz(-2.1148966) q[1];
sx q[1];
rz(-1.0245163) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56652503) q[0];
sx q[0];
rz(-1.1199513) q[0];
sx q[0];
rz(-1.2498115) q[0];
rz(-pi) q[1];
rz(2.5249135) q[2];
sx q[2];
rz(-0.70313912) q[2];
sx q[2];
rz(0.48381915) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.4752794) q[1];
sx q[1];
rz(-0.76994714) q[1];
sx q[1];
rz(-0.85789263) q[1];
x q[2];
rz(2.5173353) q[3];
sx q[3];
rz(-2.5578824) q[3];
sx q[3];
rz(2.7492085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.7496877) q[2];
sx q[2];
rz(-1.2140032) q[2];
sx q[2];
rz(2.270703) q[2];
rz(-1.332256) q[3];
sx q[3];
rz(-0.94227666) q[3];
sx q[3];
rz(1.4107305) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9695327) q[0];
sx q[0];
rz(-1.4071858) q[0];
sx q[0];
rz(0.55091888) q[0];
rz(-0.46539601) q[1];
sx q[1];
rz(-2.4702563) q[1];
sx q[1];
rz(2.904772) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2315002) q[0];
sx q[0];
rz(-1.6084558) q[0];
sx q[0];
rz(1.4715657) q[0];
rz(-pi) q[1];
rz(-2.4004164) q[2];
sx q[2];
rz(-0.20222649) q[2];
sx q[2];
rz(1.2349595) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.19141087) q[1];
sx q[1];
rz(-1.6457874) q[1];
sx q[1];
rz(-1.8412776) q[1];
rz(-0.83644609) q[3];
sx q[3];
rz(-1.340938) q[3];
sx q[3];
rz(0.072007192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.6992496) q[2];
sx q[2];
rz(-0.68168679) q[2];
sx q[2];
rz(-0.44011763) q[2];
rz(1.0951428) q[3];
sx q[3];
rz(-1.4705855) q[3];
sx q[3];
rz(1.7949036) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31297627) q[0];
sx q[0];
rz(-1.3422817) q[0];
sx q[0];
rz(-0.029504689) q[0];
rz(2.1942031) q[1];
sx q[1];
rz(-0.82059971) q[1];
sx q[1];
rz(0.11880076) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7293852) q[0];
sx q[0];
rz(-2.8280624) q[0];
sx q[0];
rz(0.14051147) q[0];
rz(-pi) q[1];
rz(-0.88840862) q[2];
sx q[2];
rz(-0.66328045) q[2];
sx q[2];
rz(-2.1036489) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5549705) q[1];
sx q[1];
rz(-1.277031) q[1];
sx q[1];
rz(0.48420669) q[1];
x q[2];
rz(-2.0649672) q[3];
sx q[3];
rz(-1.4255187) q[3];
sx q[3];
rz(0.52500099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.7416731) q[2];
sx q[2];
rz(-0.46367773) q[2];
sx q[2];
rz(1.5681533) q[2];
rz(1.9741156) q[3];
sx q[3];
rz(-1.7220595) q[3];
sx q[3];
rz(-1.8673816) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2385999) q[0];
sx q[0];
rz(-2.3165343) q[0];
sx q[0];
rz(-2.9883244) q[0];
rz(-2.080147) q[1];
sx q[1];
rz(-2.5669284) q[1];
sx q[1];
rz(1.9877888) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.518084) q[0];
sx q[0];
rz(-2.6797446) q[0];
sx q[0];
rz(-2.5005546) q[0];
rz(2.597528) q[2];
sx q[2];
rz(-0.96368921) q[2];
sx q[2];
rz(1.1407167) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-3.1187763) q[1];
sx q[1];
rz(-0.21636848) q[1];
sx q[1];
rz(1.8323891) q[1];
rz(-pi) q[2];
x q[2];
rz(0.12113573) q[3];
sx q[3];
rz(-2.2545358) q[3];
sx q[3];
rz(2.5072806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.8447421) q[2];
sx q[2];
rz(-1.9694318) q[2];
sx q[2];
rz(1.6142169) q[2];
rz(0.99689233) q[3];
sx q[3];
rz(-1.6485873) q[3];
sx q[3];
rz(-0.87866384) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-0.24755724) q[0];
sx q[0];
rz(-1.0972801) q[0];
sx q[0];
rz(-0.19009185) q[0];
rz(-2.4709573) q[1];
sx q[1];
rz(-1.1963444) q[1];
sx q[1];
rz(1.488283) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9207536) q[0];
sx q[0];
rz(-1.6016042) q[0];
sx q[0];
rz(3.0160745) q[0];
rz(-pi) q[1];
x q[1];
rz(2.699244) q[2];
sx q[2];
rz(-1.3253691) q[2];
sx q[2];
rz(-2.6128212) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.8390159) q[1];
sx q[1];
rz(-0.95514983) q[1];
sx q[1];
rz(-0.93977309) q[1];
rz(-pi) q[2];
rz(-2.2404352) q[3];
sx q[3];
rz(-2.3271051) q[3];
sx q[3];
rz(0.27064532) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.0037447475) q[2];
sx q[2];
rz(-0.90128428) q[2];
sx q[2];
rz(2.2424662) q[2];
rz(0.51504618) q[3];
sx q[3];
rz(-2.6656272) q[3];
sx q[3];
rz(-2.9639444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0638194) q[0];
sx q[0];
rz(-1.780092) q[0];
sx q[0];
rz(-0.5583981) q[0];
rz(0.28941119) q[1];
sx q[1];
rz(-0.90129539) q[1];
sx q[1];
rz(1.7064066) q[1];
rz(-2.8171956) q[2];
sx q[2];
rz(-0.73642052) q[2];
sx q[2];
rz(-3.0036075) q[2];
rz(-1.9990986) q[3];
sx q[3];
rz(-0.29853587) q[3];
sx q[3];
rz(1.5034911) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];