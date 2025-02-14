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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.762142) q[0];
sx q[0];
rz(-2.6317959) q[0];
sx q[0];
rz(2.9708746) q[0];
rz(-pi) q[1];
rz(0.40387965) q[2];
sx q[2];
rz(-0.35940659) q[2];
sx q[2];
rz(-0.071381005) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.9567555) q[1];
sx q[1];
rz(-2.7606899) q[1];
sx q[1];
rz(2.95762) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1175577) q[3];
sx q[3];
rz(-0.49542483) q[3];
sx q[3];
rz(-1.6273496) q[3];
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
rz(-1.2827778) q[3];
sx q[3];
rz(-1.1056113) q[3];
sx q[3];
rz(-0.58853308) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.823536) q[0];
sx q[0];
rz(-2.7243491) q[0];
sx q[0];
rz(2.4349924) q[0];
rz(0.63287863) q[1];
sx q[1];
rz(-2.0426079) q[1];
sx q[1];
rz(-2.852827) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5486886) q[0];
sx q[0];
rz(-0.0026772896) q[0];
sx q[0];
rz(-2.755318) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9339004) q[2];
sx q[2];
rz(-1.4121659) q[2];
sx q[2];
rz(-0.52244782) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-3.0367248) q[1];
sx q[1];
rz(-2.782208) q[1];
sx q[1];
rz(-0.64175989) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7004174) q[3];
sx q[3];
rz(-0.81038332) q[3];
sx q[3];
rz(1.2562417) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.4660945) q[2];
sx q[2];
rz(-1.0885295) q[2];
sx q[2];
rz(-1.6509854) q[2];
rz(0.7771107) q[3];
sx q[3];
rz(-1.6831574) q[3];
sx q[3];
rz(1.3597663) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31558388) q[0];
sx q[0];
rz(-1.4588139) q[0];
sx q[0];
rz(2.7050731) q[0];
rz(-0.79477683) q[1];
sx q[1];
rz(-1.7710641) q[1];
sx q[1];
rz(-2.2025542) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3600814) q[0];
sx q[0];
rz(-2.0716801) q[0];
sx q[0];
rz(-1.8003746) q[0];
rz(-1.0056507) q[2];
sx q[2];
rz(-0.88800231) q[2];
sx q[2];
rz(-1.2393774) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.446167) q[1];
sx q[1];
rz(-1.0040196) q[1];
sx q[1];
rz(-1.1935463) q[1];
x q[2];
rz(1.2235997) q[3];
sx q[3];
rz(-0.58230647) q[3];
sx q[3];
rz(-2.3045185) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.7006435) q[2];
sx q[2];
rz(-1.0627397) q[2];
sx q[2];
rz(-2.5062594) q[2];
rz(-2.3144531) q[3];
sx q[3];
rz(-2.6878036) q[3];
sx q[3];
rz(-1.2686096) q[3];
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
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24469911) q[0];
sx q[0];
rz(-3.0785955) q[0];
sx q[0];
rz(0.52198207) q[0];
rz(-0.43102795) q[1];
sx q[1];
rz(-1.6072175) q[1];
sx q[1];
rz(-0.65188754) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70744694) q[0];
sx q[0];
rz(-1.0672671) q[0];
sx q[0];
rz(-1.8412526) q[0];
rz(0.84649936) q[2];
sx q[2];
rz(-1.1870459) q[2];
sx q[2];
rz(0.64696128) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.74002162) q[1];
sx q[1];
rz(-1.5570988) q[1];
sx q[1];
rz(1.1915951) q[1];
rz(-pi) q[2];
rz(0.68199358) q[3];
sx q[3];
rz(-2.3122283) q[3];
sx q[3];
rz(-0.030205848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7151457) q[2];
sx q[2];
rz(-1.8603674) q[2];
sx q[2];
rz(2.35671) q[2];
rz(2.6134885) q[3];
sx q[3];
rz(-1.2930861) q[3];
sx q[3];
rz(1.8538792) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2454979) q[0];
sx q[0];
rz(-2.446785) q[0];
sx q[0];
rz(-2.3520663) q[0];
rz(-2.6441669) q[1];
sx q[1];
rz(-2.1010294) q[1];
sx q[1];
rz(-1.8400037) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.042846366) q[0];
sx q[0];
rz(-1.3912956) q[0];
sx q[0];
rz(2.5371206) q[0];
rz(-pi) q[1];
rz(-1.5100432) q[2];
sx q[2];
rz(-1.2323325) q[2];
sx q[2];
rz(-2.0450704) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.19986049) q[1];
sx q[1];
rz(-2.2944415) q[1];
sx q[1];
rz(1.8017215) q[1];
x q[2];
rz(0.022707247) q[3];
sx q[3];
rz(-1.2035666) q[3];
sx q[3];
rz(2.2677088) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.5975981) q[2];
sx q[2];
rz(-1.5388637) q[2];
sx q[2];
rz(0.27628118) q[2];
rz(0.9497408) q[3];
sx q[3];
rz(-0.26849982) q[3];
sx q[3];
rz(-2.5883519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62040579) q[0];
sx q[0];
rz(-0.036245417) q[0];
sx q[0];
rz(2.1976443) q[0];
rz(-1.2983407) q[1];
sx q[1];
rz(-2.0746168) q[1];
sx q[1];
rz(0.7775158) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80332887) q[0];
sx q[0];
rz(-1.9247806) q[0];
sx q[0];
rz(-0.5838809) q[0];
x q[1];
rz(-0.28892679) q[2];
sx q[2];
rz(-0.75519604) q[2];
sx q[2];
rz(0.64695219) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(3.0447415) q[1];
sx q[1];
rz(-1.6975843) q[1];
sx q[1];
rz(-0.37146588) q[1];
x q[2];
rz(3.0367682) q[3];
sx q[3];
rz(-2.8358805) q[3];
sx q[3];
rz(1.0979528) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.5139318) q[2];
sx q[2];
rz(-1.3769423) q[2];
sx q[2];
rz(-0.39109209) q[2];
rz(-0.62134653) q[3];
sx q[3];
rz(-2.4070599) q[3];
sx q[3];
rz(-0.77596107) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
sx q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4710627) q[0];
sx q[0];
rz(-3.0369018) q[0];
sx q[0];
rz(-1.4290357) q[0];
rz(2.1757226) q[1];
sx q[1];
rz(-1.254225) q[1];
sx q[1];
rz(-0.78580725) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.98322372) q[0];
sx q[0];
rz(-2.2136723) q[0];
sx q[0];
rz(1.6352562) q[0];
x q[1];
rz(1.1675535) q[2];
sx q[2];
rz(-0.87783646) q[2];
sx q[2];
rz(-0.6374109) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.2715258) q[1];
sx q[1];
rz(-1.1114239) q[1];
sx q[1];
rz(1.4794631) q[1];
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
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.4449731) q[2];
sx q[2];
rz(-0.61903054) q[2];
sx q[2];
rz(0.49717286) q[2];
rz(-0.87388006) q[3];
sx q[3];
rz(-2.0495575) q[3];
sx q[3];
rz(1.2018275) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1931605) q[0];
sx q[0];
rz(-0.30696294) q[0];
sx q[0];
rz(-0.69751414) q[0];
rz(2.7613617) q[1];
sx q[1];
rz(-1.1153328) q[1];
sx q[1];
rz(-3.0013705) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8070458) q[0];
sx q[0];
rz(-1.7522893) q[0];
sx q[0];
rz(1.6491778) q[0];
x q[1];
rz(1.5319139) q[2];
sx q[2];
rz(-2.011353) q[2];
sx q[2];
rz(-1.3300542) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.6483874) q[1];
sx q[1];
rz(-2.3304061) q[1];
sx q[1];
rz(-2.5757575) q[1];
x q[2];
rz(0.41105627) q[3];
sx q[3];
rz(-1.3945763) q[3];
sx q[3];
rz(-2.7358219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.85072881) q[2];
sx q[2];
rz(-1.2173434) q[2];
sx q[2];
rz(-2.6503837) q[2];
rz(-0.20719191) q[3];
sx q[3];
rz(-2.5184293) q[3];
sx q[3];
rz(-1.7306958) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(2.9687013) q[0];
sx q[0];
rz(-1.4145114) q[0];
sx q[0];
rz(-0.15039314) q[0];
rz(0.70101678) q[1];
sx q[1];
rz(-1.0944518) q[1];
sx q[1];
rz(1.6640123) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9532861) q[0];
sx q[0];
rz(-2.4948848) q[0];
sx q[0];
rz(-2.9394399) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4544746) q[2];
sx q[2];
rz(-1.4106299) q[2];
sx q[2];
rz(0.016906658) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.36819283) q[1];
sx q[1];
rz(-0.18562323) q[1];
sx q[1];
rz(-2.5297861) q[1];
rz(0.83563135) q[3];
sx q[3];
rz(-1.366525) q[3];
sx q[3];
rz(-2.6428619) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.574934) q[2];
sx q[2];
rz(-2.7329972) q[2];
sx q[2];
rz(-1.6179786) q[2];
rz(-2.5426215) q[3];
sx q[3];
rz(-0.50061575) q[3];
sx q[3];
rz(-0.18794255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4481675) q[0];
sx q[0];
rz(-0.42891112) q[0];
sx q[0];
rz(-1.6397788) q[0];
rz(2.2046454) q[1];
sx q[1];
rz(-1.0277156) q[1];
sx q[1];
rz(-2.1243336) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1382683) q[0];
sx q[0];
rz(-1.956245) q[0];
sx q[0];
rz(1.9786563) q[0];
rz(2.4577519) q[2];
sx q[2];
rz(-0.13540395) q[2];
sx q[2];
rz(1.1410448) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.63314519) q[1];
sx q[1];
rz(-0.26403759) q[1];
sx q[1];
rz(1.6980843) q[1];
rz(0.42041619) q[3];
sx q[3];
rz(-0.46483332) q[3];
sx q[3];
rz(0.55373389) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0246058) q[2];
sx q[2];
rz(-2.459343) q[2];
sx q[2];
rz(1.619722) q[2];
rz(-1.4874602) q[3];
sx q[3];
rz(-1.6493075) q[3];
sx q[3];
rz(1.7447932) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7617154) q[0];
sx q[0];
rz(-1.3991671) q[0];
sx q[0];
rz(0.73200926) q[0];
rz(-1.6666182) q[1];
sx q[1];
rz(-1.9261618) q[1];
sx q[1];
rz(0.7934657) q[1];
rz(-0.39948612) q[2];
sx q[2];
rz(-2.1151092) q[2];
sx q[2];
rz(-1.5110037) q[2];
rz(2.8216437) q[3];
sx q[3];
rz(-1.2816164) q[3];
sx q[3];
rz(1.5641227) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
