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
rz(0.30492914) q[0];
sx q[0];
rz(-0.066216901) q[0];
sx q[0];
rz(-2.9856292) q[0];
rz(1.1671542) q[1];
sx q[1];
rz(3.7937556) q[1];
sx q[1];
rz(8.9383386) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64336813) q[0];
sx q[0];
rz(-2.273145) q[0];
sx q[0];
rz(0.81612103) q[0];
rz(-0.90801348) q[2];
sx q[2];
rz(-2.5129299) q[2];
sx q[2];
rz(-2.6952621) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.26775441) q[1];
sx q[1];
rz(-0.43711409) q[1];
sx q[1];
rz(0.74437352) q[1];
rz(-pi) q[2];
rz(-0.82709978) q[3];
sx q[3];
rz(-1.4129644) q[3];
sx q[3];
rz(-2.1438897) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.277694) q[2];
sx q[2];
rz(-2.4319477) q[2];
sx q[2];
rz(2.7126183) q[2];
rz(-0.046836827) q[3];
sx q[3];
rz(-2.7497079) q[3];
sx q[3];
rz(-0.61280167) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.866339) q[0];
sx q[0];
rz(-2.1493122) q[0];
sx q[0];
rz(0.66542768) q[0];
rz(1.1412507) q[1];
sx q[1];
rz(-1.7612061) q[1];
sx q[1];
rz(0.64249396) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2739853) q[0];
sx q[0];
rz(-1.6216177) q[0];
sx q[0];
rz(-2.4073868) q[0];
x q[1];
rz(0.69509721) q[2];
sx q[2];
rz(-1.8029658) q[2];
sx q[2];
rz(2.1921981) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.2765668) q[1];
sx q[1];
rz(-1.5986018) q[1];
sx q[1];
rz(0.69827484) q[1];
x q[2];
rz(-0.1763686) q[3];
sx q[3];
rz(-2.1049728) q[3];
sx q[3];
rz(1.4458223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.3147754) q[2];
sx q[2];
rz(-2.3544957) q[2];
sx q[2];
rz(-0.55244201) q[2];
rz(1.6636482) q[3];
sx q[3];
rz(-1.843957) q[3];
sx q[3];
rz(-2.4655931) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8058572) q[0];
sx q[0];
rz(-2.4215846) q[0];
sx q[0];
rz(2.3190401) q[0];
rz(-2.3303253) q[1];
sx q[1];
rz(-0.30458105) q[1];
sx q[1];
rz(1.4623581) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84483355) q[0];
sx q[0];
rz(-2.5243596) q[0];
sx q[0];
rz(-2.4485641) q[0];
x q[1];
rz(0.16872318) q[2];
sx q[2];
rz(-2.7232735) q[2];
sx q[2];
rz(0.63101286) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9595527) q[1];
sx q[1];
rz(-2.5904852) q[1];
sx q[1];
rz(-2.1564547) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.30841804) q[3];
sx q[3];
rz(-1.8659411) q[3];
sx q[3];
rz(0.67527387) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2751969) q[2];
sx q[2];
rz(-0.82922816) q[2];
sx q[2];
rz(-2.5466476) q[2];
rz(-2.7368937) q[3];
sx q[3];
rz(-2.0345104) q[3];
sx q[3];
rz(1.2987202) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73612708) q[0];
sx q[0];
rz(-1.8647702) q[0];
sx q[0];
rz(-0.24293105) q[0];
rz(-2.2566707) q[1];
sx q[1];
rz(-0.72598571) q[1];
sx q[1];
rz(0.0028217908) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.045319917) q[0];
sx q[0];
rz(-2.3191873) q[0];
sx q[0];
rz(-2.1069645) q[0];
x q[1];
rz(-2.1921333) q[2];
sx q[2];
rz(-2.0120125) q[2];
sx q[2];
rz(-0.37829933) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1403919) q[1];
sx q[1];
rz(-2.2691548) q[1];
sx q[1];
rz(-3.0496518) q[1];
x q[2];
rz(-0.61496998) q[3];
sx q[3];
rz(-2.0738064) q[3];
sx q[3];
rz(0.78395432) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.60488492) q[2];
sx q[2];
rz(-1.1515836) q[2];
sx q[2];
rz(2.4082157) q[2];
rz(0.65507656) q[3];
sx q[3];
rz(-0.30787444) q[3];
sx q[3];
rz(-0.56183279) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4578399) q[0];
sx q[0];
rz(-0.34128749) q[0];
sx q[0];
rz(-0.14232464) q[0];
rz(-1.3617474) q[1];
sx q[1];
rz(-0.35265499) q[1];
sx q[1];
rz(0.49837643) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67923421) q[0];
sx q[0];
rz(-0.7158044) q[0];
sx q[0];
rz(1.6648034) q[0];
rz(-pi) q[1];
x q[1];
rz(0.054008897) q[2];
sx q[2];
rz(-2.3819792) q[2];
sx q[2];
rz(0.97689607) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.315334) q[1];
sx q[1];
rz(-2.3439601) q[1];
sx q[1];
rz(-2.7263097) q[1];
rz(0.33739319) q[3];
sx q[3];
rz(-2.3327347) q[3];
sx q[3];
rz(-1.0800501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0087697) q[2];
sx q[2];
rz(-1.885773) q[2];
sx q[2];
rz(-0.69457561) q[2];
rz(-0.83235598) q[3];
sx q[3];
rz(-1.1135626) q[3];
sx q[3];
rz(0.25025234) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4271127) q[0];
sx q[0];
rz(-2.6321593) q[0];
sx q[0];
rz(-0.66202778) q[0];
rz(-1.1854677) q[1];
sx q[1];
rz(-1.0292425) q[1];
sx q[1];
rz(-1.9756636) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53138779) q[0];
sx q[0];
rz(-2.234196) q[0];
sx q[0];
rz(-1.8276455) q[0];
x q[1];
rz(-1.7244965) q[2];
sx q[2];
rz(-0.96818334) q[2];
sx q[2];
rz(-0.034687925) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.1665713) q[1];
sx q[1];
rz(-1.6432883) q[1];
sx q[1];
rz(0.049866421) q[1];
x q[2];
rz(-0.95997973) q[3];
sx q[3];
rz(-1.0418384) q[3];
sx q[3];
rz(2.4774266) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.7621496) q[2];
sx q[2];
rz(-1.9337485) q[2];
sx q[2];
rz(0.70972788) q[2];
rz(-0.48192561) q[3];
sx q[3];
rz(-2.66633) q[3];
sx q[3];
rz(-0.019439241) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67126453) q[0];
sx q[0];
rz(-2.1605991) q[0];
sx q[0];
rz(-1.5671267) q[0];
rz(1.6471242) q[1];
sx q[1];
rz(-0.44691214) q[1];
sx q[1];
rz(-0.87237298) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39473662) q[0];
sx q[0];
rz(-2.7867705) q[0];
sx q[0];
rz(1.0635934) q[0];
x q[1];
rz(-2.4192737) q[2];
sx q[2];
rz(-1.7120984) q[2];
sx q[2];
rz(1.5324805) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.38546121) q[1];
sx q[1];
rz(-1.4646936) q[1];
sx q[1];
rz(-0.030563449) q[1];
rz(3.116901) q[3];
sx q[3];
rz(-2.2253312) q[3];
sx q[3];
rz(0.86658044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.81277728) q[2];
sx q[2];
rz(-0.8605364) q[2];
sx q[2];
rz(-1.1278197) q[2];
rz(-0.40337107) q[3];
sx q[3];
rz(-1.4453459) q[3];
sx q[3];
rz(2.4283714) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.23680747) q[0];
sx q[0];
rz(-1.1431563) q[0];
sx q[0];
rz(-1.8444201) q[0];
rz(2.6203268) q[1];
sx q[1];
rz(-1.2626941) q[1];
sx q[1];
rz(-0.062779471) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.990898) q[0];
sx q[0];
rz(-0.44838649) q[0];
sx q[0];
rz(1.4175116) q[0];
rz(-pi) q[1];
rz(-0.36022236) q[2];
sx q[2];
rz(-1.3531605) q[2];
sx q[2];
rz(1.9955903) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.081689) q[1];
sx q[1];
rz(-1.4270596) q[1];
sx q[1];
rz(1.8200851) q[1];
rz(-pi) q[2];
x q[2];
rz(2.217123) q[3];
sx q[3];
rz(-1.1271584) q[3];
sx q[3];
rz(-2.8868589) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.2433743) q[2];
sx q[2];
rz(-2.2369907) q[2];
sx q[2];
rz(2.173219) q[2];
rz(2.6280256) q[3];
sx q[3];
rz(-1.3526724) q[3];
sx q[3];
rz(0.33094049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53720713) q[0];
sx q[0];
rz(-2.4810915) q[0];
sx q[0];
rz(-0.78395098) q[0];
rz(-1.9369269) q[1];
sx q[1];
rz(-0.37310633) q[1];
sx q[1];
rz(-1.3752259) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.169302) q[0];
sx q[0];
rz(-1.8344886) q[0];
sx q[0];
rz(1.3486116) q[0];
rz(-3.1143931) q[2];
sx q[2];
rz(-1.4157032) q[2];
sx q[2];
rz(-1.9981249) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.7837401) q[1];
sx q[1];
rz(-2.2688818) q[1];
sx q[1];
rz(0.74905209) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9638644) q[3];
sx q[3];
rz(-1.8788218) q[3];
sx q[3];
rz(0.34402381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.2969926) q[2];
sx q[2];
rz(-2.832909) q[2];
sx q[2];
rz(2.6958579) q[2];
rz(-0.60278696) q[3];
sx q[3];
rz(-1.2647537) q[3];
sx q[3];
rz(0.84356892) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26850253) q[0];
sx q[0];
rz(-0.52274811) q[0];
sx q[0];
rz(-2.6173746) q[0];
rz(1.9408608) q[1];
sx q[1];
rz(-1.8914765) q[1];
sx q[1];
rz(-1.1842309) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6551334) q[0];
sx q[0];
rz(-2.7441027) q[0];
sx q[0];
rz(-2.6349806) q[0];
rz(-pi) q[1];
rz(-1.7767364) q[2];
sx q[2];
rz(-2.8194234) q[2];
sx q[2];
rz(1.8274771) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.56630707) q[1];
sx q[1];
rz(-0.13338365) q[1];
sx q[1];
rz(-3.0493983) q[1];
rz(-pi) q[2];
rz(-3.0940787) q[3];
sx q[3];
rz(-2.5121452) q[3];
sx q[3];
rz(0.87241064) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.9253917) q[2];
sx q[2];
rz(-0.4500469) q[2];
sx q[2];
rz(2.0897384) q[2];
rz(-0.22107407) q[3];
sx q[3];
rz(-1.3007921) q[3];
sx q[3];
rz(2.2033447) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5939519) q[0];
sx q[0];
rz(-1.6407536) q[0];
sx q[0];
rz(-3.092691) q[0];
rz(-0.37638695) q[1];
sx q[1];
rz(-2.2793437) q[1];
sx q[1];
rz(1.9752165) q[1];
rz(0.049565774) q[2];
sx q[2];
rz(-0.57065806) q[2];
sx q[2];
rz(-1.5869637) q[2];
rz(2.0509649) q[3];
sx q[3];
rz(-2.2259897) q[3];
sx q[3];
rz(2.5799354) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
