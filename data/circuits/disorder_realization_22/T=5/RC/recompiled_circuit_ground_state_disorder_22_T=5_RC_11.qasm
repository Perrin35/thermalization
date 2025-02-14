OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.8393505) q[0];
sx q[0];
rz(-1.8678764) q[0];
sx q[0];
rz(-0.94925517) q[0];
rz(-1.0957837) q[1];
sx q[1];
rz(2.1646808) q[1];
sx q[1];
rz(10.77471) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5744517) q[0];
sx q[0];
rz(-2.0724839) q[0];
sx q[0];
rz(-1.4760963) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4241997) q[2];
sx q[2];
rz(-1.9001385) q[2];
sx q[2];
rz(2.6419287) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.9266859) q[1];
sx q[1];
rz(-1.6388571) q[1];
sx q[1];
rz(-2.7665274) q[1];
x q[2];
rz(-2.0032521) q[3];
sx q[3];
rz(-1.3210332) q[3];
sx q[3];
rz(-2.7065211) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.533941) q[2];
sx q[2];
rz(-2.0484296) q[2];
sx q[2];
rz(-0.1926113) q[2];
rz(1.8588148) q[3];
sx q[3];
rz(-1.1056113) q[3];
sx q[3];
rz(-0.58853308) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31805661) q[0];
sx q[0];
rz(-2.7243491) q[0];
sx q[0];
rz(0.70660025) q[0];
rz(-2.508714) q[1];
sx q[1];
rz(-2.0426079) q[1];
sx q[1];
rz(-2.852827) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59161883) q[0];
sx q[0];
rz(-1.571805) q[0];
sx q[0];
rz(0.002480025) q[0];
x q[1];
rz(-2.9339004) q[2];
sx q[2];
rz(-1.7294267) q[2];
sx q[2];
rz(0.52244782) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.10486785) q[1];
sx q[1];
rz(-0.35938469) q[1];
sx q[1];
rz(2.4998328) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.4411753) q[3];
sx q[3];
rz(-2.3312093) q[3];
sx q[3];
rz(1.885351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.67549813) q[2];
sx q[2];
rz(-1.0885295) q[2];
sx q[2];
rz(-1.6509854) q[2];
rz(2.364482) q[3];
sx q[3];
rz(-1.6831574) q[3];
sx q[3];
rz(-1.3597663) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31558388) q[0];
sx q[0];
rz(-1.4588139) q[0];
sx q[0];
rz(2.7050731) q[0];
rz(0.79477683) q[1];
sx q[1];
rz(-1.7710641) q[1];
sx q[1];
rz(-0.93903843) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3600814) q[0];
sx q[0];
rz(-1.0699125) q[0];
sx q[0];
rz(1.341218) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.58231488) q[2];
sx q[2];
rz(-2.2852201) q[2];
sx q[2];
rz(-1.1143052) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.3311829) q[1];
sx q[1];
rz(-0.66920921) q[1];
sx q[1];
rz(-2.616908) q[1];
rz(-pi) q[2];
x q[2];
rz(2.125199) q[3];
sx q[3];
rz(-1.3825584) q[3];
sx q[3];
rz(0.44023006) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.44094917) q[2];
sx q[2];
rz(-2.078853) q[2];
sx q[2];
rz(2.5062594) q[2];
rz(-0.82713953) q[3];
sx q[3];
rz(-0.45378903) q[3];
sx q[3];
rz(-1.2686096) q[3];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8968935) q[0];
sx q[0];
rz(-0.06299717) q[0];
sx q[0];
rz(-0.52198207) q[0];
rz(-0.43102795) q[1];
sx q[1];
rz(-1.5343752) q[1];
sx q[1];
rz(0.65188754) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7303607) q[0];
sx q[0];
rz(-1.8069977) q[0];
sx q[0];
rz(-0.51933164) q[0];
rz(-pi) q[1];
rz(-2.6471557) q[2];
sx q[2];
rz(-2.2324413) q[2];
sx q[2];
rz(-2.5376157) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.3451772) q[1];
sx q[1];
rz(-2.762156) q[1];
sx q[1];
rz(-1.5338083) q[1];
x q[2];
rz(0.96792696) q[3];
sx q[3];
rz(-0.9612007) q[3];
sx q[3];
rz(2.2948007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.7151457) q[2];
sx q[2];
rz(-1.2812252) q[2];
sx q[2];
rz(-0.78488266) q[2];
rz(0.52810413) q[3];
sx q[3];
rz(-1.2930861) q[3];
sx q[3];
rz(-1.8538792) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89609471) q[0];
sx q[0];
rz(-2.446785) q[0];
sx q[0];
rz(-2.3520663) q[0];
rz(0.49742571) q[1];
sx q[1];
rz(-2.1010294) q[1];
sx q[1];
rz(-1.8400037) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7808944) q[0];
sx q[0];
rz(-0.62735617) q[0];
sx q[0];
rz(-0.30904667) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6315494) q[2];
sx q[2];
rz(-1.9092602) q[2];
sx q[2];
rz(-2.0450704) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.19986049) q[1];
sx q[1];
rz(-0.84715119) q[1];
sx q[1];
rz(1.8017215) q[1];
rz(1.6297518) q[3];
sx q[3];
rz(-0.36789933) q[3];
sx q[3];
rz(-2.2045362) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.54399458) q[2];
sx q[2];
rz(-1.5388637) q[2];
sx q[2];
rz(0.27628118) q[2];
rz(-0.9497408) q[3];
sx q[3];
rz(-0.26849982) q[3];
sx q[3];
rz(2.5883519) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5211869) q[0];
sx q[0];
rz(-0.036245417) q[0];
sx q[0];
rz(0.94394839) q[0];
rz(1.2983407) q[1];
sx q[1];
rz(-1.0669758) q[1];
sx q[1];
rz(-2.3640769) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28431129) q[0];
sx q[0];
rz(-2.4696283) q[0];
sx q[0];
rz(-2.551033) q[0];
rz(-pi) q[1];
rz(-0.73410122) q[2];
sx q[2];
rz(-1.3742374) q[2];
sx q[2];
rz(2.4308956) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.4247269) q[1];
sx q[1];
rz(-1.9391372) q[1];
sx q[1];
rz(-1.7067545) q[1];
rz(1.603807) q[3];
sx q[3];
rz(-1.8747765) q[3];
sx q[3];
rz(-2.153521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.5139318) q[2];
sx q[2];
rz(-1.3769423) q[2];
sx q[2];
rz(0.39109209) q[2];
rz(-0.62134653) q[3];
sx q[3];
rz(-0.73453271) q[3];
sx q[3];
rz(-2.3656316) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4710627) q[0];
sx q[0];
rz(-0.10469086) q[0];
sx q[0];
rz(1.4290357) q[0];
rz(2.1757226) q[1];
sx q[1];
rz(-1.8873676) q[1];
sx q[1];
rz(-2.3557854) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62625058) q[0];
sx q[0];
rz(-1.6223755) q[0];
sx q[0];
rz(-0.64387384) q[0];
rz(-pi) q[1];
x q[1];
rz(0.73432095) q[2];
sx q[2];
rz(-1.2641126) q[2];
sx q[2];
rz(-1.9421362) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.475226) q[1];
sx q[1];
rz(-0.46772784) q[1];
sx q[1];
rz(-0.18233129) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8604516) q[3];
sx q[3];
rz(-1.9928586) q[3];
sx q[3];
rz(0.45988032) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4449731) q[2];
sx q[2];
rz(-2.5225621) q[2];
sx q[2];
rz(2.6444198) q[2];
rz(2.2677126) q[3];
sx q[3];
rz(-2.0495575) q[3];
sx q[3];
rz(1.2018275) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1931605) q[0];
sx q[0];
rz(-2.8346297) q[0];
sx q[0];
rz(-0.69751414) q[0];
rz(-0.3802309) q[1];
sx q[1];
rz(-2.0262599) q[1];
sx q[1];
rz(-0.14022216) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8070458) q[0];
sx q[0];
rz(-1.7522893) q[0];
sx q[0];
rz(-1.6491778) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.082265286) q[2];
sx q[2];
rz(-2.6994355) q[2];
sx q[2];
rz(-1.7205659) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6517087) q[1];
sx q[1];
rz(-1.1715285) q[1];
sx q[1];
rz(-0.7266161) q[1];
x q[2];
rz(-1.7626552) q[3];
sx q[3];
rz(-1.975112) q[3];
sx q[3];
rz(-1.9003001) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.2908638) q[2];
sx q[2];
rz(-1.9242492) q[2];
sx q[2];
rz(2.6503837) q[2];
rz(0.20719191) q[3];
sx q[3];
rz(-0.62316337) q[3];
sx q[3];
rz(1.4108968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9687013) q[0];
sx q[0];
rz(-1.7270813) q[0];
sx q[0];
rz(2.9911995) q[0];
rz(-2.4405759) q[1];
sx q[1];
rz(-2.0471408) q[1];
sx q[1];
rz(1.4775803) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9532861) q[0];
sx q[0];
rz(-0.64670783) q[0];
sx q[0];
rz(0.20215277) q[0];
rz(-pi) q[1];
rz(-1.4544746) q[2];
sx q[2];
rz(-1.4106299) q[2];
sx q[2];
rz(-0.016906658) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.3353053) q[1];
sx q[1];
rz(-1.6769969) q[1];
sx q[1];
rz(0.152529) q[1];
rz(-pi) q[2];
rz(-1.870369) q[3];
sx q[3];
rz(-2.3837187) q[3];
sx q[3];
rz(-1.8488499) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.56665862) q[2];
sx q[2];
rz(-2.7329972) q[2];
sx q[2];
rz(1.5236141) q[2];
rz(-0.59897113) q[3];
sx q[3];
rz(-2.6409769) q[3];
sx q[3];
rz(-0.18794255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69342518) q[0];
sx q[0];
rz(-0.42891112) q[0];
sx q[0];
rz(-1.6397788) q[0];
rz(-0.93694726) q[1];
sx q[1];
rz(-2.1138771) q[1];
sx q[1];
rz(-1.017259) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0033244) q[0];
sx q[0];
rz(-1.956245) q[0];
sx q[0];
rz(1.1629364) q[0];
rz(-pi) q[1];
rz(0.10521484) q[2];
sx q[2];
rz(-1.485409) q[2];
sx q[2];
rz(2.8920023) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.76495095) q[1];
sx q[1];
rz(-1.8326474) q[1];
sx q[1];
rz(0.034305926) q[1];
rz(-pi) q[2];
rz(-2.7211765) q[3];
sx q[3];
rz(-2.6767593) q[3];
sx q[3];
rz(-0.55373389) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-3.0246058) q[2];
sx q[2];
rz(-2.459343) q[2];
sx q[2];
rz(-1.5218706) q[2];
rz(-1.6541325) q[3];
sx q[3];
rz(-1.6493075) q[3];
sx q[3];
rz(1.3967995) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
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
rz(2.7617154) q[0];
sx q[0];
rz(-1.3991671) q[0];
sx q[0];
rz(0.73200926) q[0];
rz(1.6666182) q[1];
sx q[1];
rz(-1.2154308) q[1];
sx q[1];
rz(-2.348127) q[1];
rz(0.99967069) q[2];
sx q[2];
rz(-0.66304211) q[2];
sx q[2];
rz(-2.1950051) q[2];
rz(2.3839604) q[3];
sx q[3];
rz(-0.42790596) q[3];
sx q[3];
rz(-0.71732646) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
