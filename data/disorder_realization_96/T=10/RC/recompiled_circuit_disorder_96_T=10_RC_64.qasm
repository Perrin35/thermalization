OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-3.0918026) q[0];
sx q[0];
rz(-3.0135305) q[0];
sx q[0];
rz(-0.81737104) q[0];
rz(0.983239) q[1];
sx q[1];
rz(-0.53951889) q[1];
sx q[1];
rz(-1.2004381) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.565633) q[0];
sx q[0];
rz(-2.6607249) q[0];
sx q[0];
rz(-2.5984882) q[0];
rz(-pi) q[1];
rz(3.009216) q[2];
sx q[2];
rz(-1.0356324) q[2];
sx q[2];
rz(-1.7420499) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.8564954) q[1];
sx q[1];
rz(-2.3772117) q[1];
sx q[1];
rz(2.7544423) q[1];
rz(0.31906268) q[3];
sx q[3];
rz(-2.3205119) q[3];
sx q[3];
rz(1.0533489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.0212705) q[2];
sx q[2];
rz(-2.5734084) q[2];
sx q[2];
rz(-1.583064) q[2];
rz(2.1448686) q[3];
sx q[3];
rz(-0.45209) q[3];
sx q[3];
rz(2.7157917) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5581756) q[0];
sx q[0];
rz(-0.75220627) q[0];
sx q[0];
rz(-3.0875207) q[0];
rz(1.1955098) q[1];
sx q[1];
rz(-1.0369438) q[1];
sx q[1];
rz(0.53584677) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5862522) q[0];
sx q[0];
rz(-0.99389168) q[0];
sx q[0];
rz(-0.14848407) q[0];
rz(-pi) q[1];
rz(2.8500154) q[2];
sx q[2];
rz(-2.4107286) q[2];
sx q[2];
rz(-0.72327327) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.6260813) q[1];
sx q[1];
rz(-1.9768081) q[1];
sx q[1];
rz(1.8477693) q[1];
rz(-pi) q[2];
rz(0.032580839) q[3];
sx q[3];
rz(-1.4279281) q[3];
sx q[3];
rz(0.25861614) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.0216996) q[2];
sx q[2];
rz(-1.0571486) q[2];
sx q[2];
rz(-2.1832441) q[2];
rz(-0.066453233) q[3];
sx q[3];
rz(-1.5840014) q[3];
sx q[3];
rz(-0.44979969) q[3];
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
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7217343) q[0];
sx q[0];
rz(-1.2092713) q[0];
sx q[0];
rz(0.15047519) q[0];
rz(-0.45723215) q[1];
sx q[1];
rz(-0.91232863) q[1];
sx q[1];
rz(-0.025807468) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8919825) q[0];
sx q[0];
rz(-1.5907856) q[0];
sx q[0];
rz(1.32094) q[0];
rz(-pi) q[1];
x q[1];
rz(1.259272) q[2];
sx q[2];
rz(-1.6561964) q[2];
sx q[2];
rz(2.4740919) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.4986213) q[1];
sx q[1];
rz(-1.1585304) q[1];
sx q[1];
rz(0.69795124) q[1];
rz(-pi) q[2];
rz(-1.8886391) q[3];
sx q[3];
rz(-0.40099537) q[3];
sx q[3];
rz(-0.81438118) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1228483) q[2];
sx q[2];
rz(-0.3652502) q[2];
sx q[2];
rz(0.5853816) q[2];
rz(-0.18150005) q[3];
sx q[3];
rz(-1.8094962) q[3];
sx q[3];
rz(-1.5649149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.240775) q[0];
sx q[0];
rz(-0.62269354) q[0];
sx q[0];
rz(0.17661072) q[0];
rz(2.2606842) q[1];
sx q[1];
rz(-2.0753588) q[1];
sx q[1];
rz(-2.6054629) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25585184) q[0];
sx q[0];
rz(-0.8937853) q[0];
sx q[0];
rz(0.52744249) q[0];
x q[1];
rz(-0.23668004) q[2];
sx q[2];
rz(-1.6263905) q[2];
sx q[2];
rz(0.4586763) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.22524658) q[1];
sx q[1];
rz(-2.1919474) q[1];
sx q[1];
rz(1.1017373) q[1];
rz(-pi) q[2];
rz(-2.7809308) q[3];
sx q[3];
rz(-2.4324527) q[3];
sx q[3];
rz(-0.053645596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.6716016) q[2];
sx q[2];
rz(-1.4208379) q[2];
sx q[2];
rz(-2.0969351) q[2];
rz(0.70703834) q[3];
sx q[3];
rz(-0.98393327) q[3];
sx q[3];
rz(2.7846591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2351284) q[0];
sx q[0];
rz(-1.3236073) q[0];
sx q[0];
rz(-1.0513603) q[0];
rz(-1.6479187) q[1];
sx q[1];
rz(-2.5525679) q[1];
sx q[1];
rz(0.043118127) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2252786) q[0];
sx q[0];
rz(-1.3818704) q[0];
sx q[0];
rz(0.96154763) q[0];
x q[1];
rz(1.0518603) q[2];
sx q[2];
rz(-1.3216002) q[2];
sx q[2];
rz(0.11617004) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.9777898) q[1];
sx q[1];
rz(-1.8352574) q[1];
sx q[1];
rz(-2.8603641) q[1];
x q[2];
rz(-2.352036) q[3];
sx q[3];
rz(-1.6766747) q[3];
sx q[3];
rz(1.7136128) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.12895) q[2];
sx q[2];
rz(-0.49762112) q[2];
sx q[2];
rz(0.35219231) q[2];
rz(-0.59018618) q[3];
sx q[3];
rz(-2.6679109) q[3];
sx q[3];
rz(0.56110704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6234289) q[0];
sx q[0];
rz(-1.9264899) q[0];
sx q[0];
rz(1.9859001) q[0];
rz(-2.39134) q[1];
sx q[1];
rz(-2.2022088) q[1];
sx q[1];
rz(2.0828784) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1832702) q[0];
sx q[0];
rz(-2.000862) q[0];
sx q[0];
rz(-1.8002611) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2054339) q[2];
sx q[2];
rz(-2.0632671) q[2];
sx q[2];
rz(-1.1577215) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.2211654) q[1];
sx q[1];
rz(-1.3472918) q[1];
sx q[1];
rz(-1.2616874) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3771463) q[3];
sx q[3];
rz(-2.2691257) q[3];
sx q[3];
rz(0.26102548) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.47026149) q[2];
sx q[2];
rz(-1.7241717) q[2];
sx q[2];
rz(-1.2188101) q[2];
rz(-1.9865215) q[3];
sx q[3];
rz(-0.23854908) q[3];
sx q[3];
rz(-1.7003805) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.041615151) q[0];
sx q[0];
rz(-1.8247373) q[0];
sx q[0];
rz(1.6301427) q[0];
rz(1.7639683) q[1];
sx q[1];
rz(-0.310985) q[1];
sx q[1];
rz(-0.84164936) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8182939) q[0];
sx q[0];
rz(-0.57384402) q[0];
sx q[0];
rz(-2.7159575) q[0];
rz(-pi) q[1];
rz(0.62692554) q[2];
sx q[2];
rz(-1.4142087) q[2];
sx q[2];
rz(0.92948929) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.6262498) q[1];
sx q[1];
rz(-1.6193577) q[1];
sx q[1];
rz(-1.4944782) q[1];
x q[2];
rz(-1.3326725) q[3];
sx q[3];
rz(-0.75546414) q[3];
sx q[3];
rz(1.8734224) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.001361751) q[2];
sx q[2];
rz(-1.3683687) q[2];
sx q[2];
rz(-0.049953071) q[2];
rz(-2.4800381) q[3];
sx q[3];
rz(-2.619132) q[3];
sx q[3];
rz(-2.9522827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2426303) q[0];
sx q[0];
rz(-1.0268509) q[0];
sx q[0];
rz(1.4021953) q[0];
rz(-0.095480355) q[1];
sx q[1];
rz(-1.9752558) q[1];
sx q[1];
rz(-0.41762525) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2554889) q[0];
sx q[0];
rz(-2.4574453) q[0];
sx q[0];
rz(0.72649254) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.278644) q[2];
sx q[2];
rz(-0.38049618) q[2];
sx q[2];
rz(0.7064864) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.8199181) q[1];
sx q[1];
rz(-0.770861) q[1];
sx q[1];
rz(0.8701156) q[1];
rz(-pi) q[2];
rz(1.1975708) q[3];
sx q[3];
rz(-1.8456568) q[3];
sx q[3];
rz(-1.1653792) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3461356) q[2];
sx q[2];
rz(-1.6980349) q[2];
sx q[2];
rz(2.0987089) q[2];
rz(-0.67388326) q[3];
sx q[3];
rz(-1.6582812) q[3];
sx q[3];
rz(-0.90014443) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3089356) q[0];
sx q[0];
rz(-1.4082264) q[0];
sx q[0];
rz(2.5119264) q[0];
rz(-2.5667403) q[1];
sx q[1];
rz(-1.31153) q[1];
sx q[1];
rz(0.94690698) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1614721) q[0];
sx q[0];
rz(-1.8959909) q[0];
sx q[0];
rz(-1.253771) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7526555) q[2];
sx q[2];
rz(-3.0367594) q[2];
sx q[2];
rz(-0.36738415) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.2448333) q[1];
sx q[1];
rz(-1.6549126) q[1];
sx q[1];
rz(-0.33478488) q[1];
rz(-pi) q[2];
rz(-0.19091786) q[3];
sx q[3];
rz(-1.1357765) q[3];
sx q[3];
rz(0.85489475) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.5809014) q[2];
sx q[2];
rz(-2.0118482) q[2];
sx q[2];
rz(-1.2488731) q[2];
rz(-2.4272264) q[3];
sx q[3];
rz(-1.2759821) q[3];
sx q[3];
rz(-2.8760288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
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
rz(-2.0666075) q[0];
sx q[0];
rz(-2.5191436) q[0];
sx q[0];
rz(2.4865436) q[0];
rz(2.24522) q[1];
sx q[1];
rz(-2.2348149) q[1];
sx q[1];
rz(2.4972829) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6039186) q[0];
sx q[0];
rz(-2.9593421) q[0];
sx q[0];
rz(2.8331579) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2868025) q[2];
sx q[2];
rz(-1.7456747) q[2];
sx q[2];
rz(1.0027494) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.555968) q[1];
sx q[1];
rz(-1.854419) q[1];
sx q[1];
rz(-0.3507627) q[1];
rz(-0.44390042) q[3];
sx q[3];
rz(-2.0022941) q[3];
sx q[3];
rz(2.3865226) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.3506938) q[2];
sx q[2];
rz(-1.6692946) q[2];
sx q[2];
rz(-2.541686) q[2];
rz(2.24263) q[3];
sx q[3];
rz(-0.18342429) q[3];
sx q[3];
rz(1.3658587) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8469289) q[0];
sx q[0];
rz(-1.3273205) q[0];
sx q[0];
rz(2.474665) q[0];
rz(0.22944336) q[1];
sx q[1];
rz(-0.89090092) q[1];
sx q[1];
rz(0.13577239) q[1];
rz(1.3953801) q[2];
sx q[2];
rz(-0.63175628) q[2];
sx q[2];
rz(0.14887688) q[2];
rz(2.7502144) q[3];
sx q[3];
rz(-2.2073675) q[3];
sx q[3];
rz(-1.7195306) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
