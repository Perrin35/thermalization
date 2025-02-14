OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.71830463) q[0];
sx q[0];
rz(-0.77646065) q[0];
sx q[0];
rz(-0.36592308) q[0];
rz(0.34691063) q[1];
sx q[1];
rz(-0.28471714) q[1];
sx q[1];
rz(1.7528344) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58908868) q[0];
sx q[0];
rz(-2.9285746) q[0];
sx q[0];
rz(0.32040839) q[0];
rz(-0.35438673) q[2];
sx q[2];
rz(-0.72438188) q[2];
sx q[2];
rz(-0.61930932) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.5386913) q[1];
sx q[1];
rz(-0.60990342) q[1];
sx q[1];
rz(-2.2944488) q[1];
x q[2];
rz(-2.5479814) q[3];
sx q[3];
rz(-2.0214012) q[3];
sx q[3];
rz(3.1385147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.1846788) q[2];
sx q[2];
rz(-0.92014402) q[2];
sx q[2];
rz(-0.93519768) q[2];
rz(-2.8397371) q[3];
sx q[3];
rz(-1.5993092) q[3];
sx q[3];
rz(-0.39366084) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97717706) q[0];
sx q[0];
rz(-1.875832) q[0];
sx q[0];
rz(-0.49829495) q[0];
rz(-1.4277108) q[1];
sx q[1];
rz(-1.9352813) q[1];
sx q[1];
rz(-1.0867585) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3918864) q[0];
sx q[0];
rz(-1.8072309) q[0];
sx q[0];
rz(2.4357585) q[0];
x q[1];
rz(-0.94765122) q[2];
sx q[2];
rz(-0.32440475) q[2];
sx q[2];
rz(2.0632921) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.9683428) q[1];
sx q[1];
rz(-2.0054711) q[1];
sx q[1];
rz(-1.3938851) q[1];
x q[2];
rz(1.558775) q[3];
sx q[3];
rz(-1.88899) q[3];
sx q[3];
rz(-2.1592922) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.19865092) q[2];
sx q[2];
rz(-0.31060878) q[2];
sx q[2];
rz(-1.8093713) q[2];
rz(1.9940935) q[3];
sx q[3];
rz(-1.56286) q[3];
sx q[3];
rz(1.8089627) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9320817) q[0];
sx q[0];
rz(-1.7387583) q[0];
sx q[0];
rz(0.15723666) q[0];
rz(-1.9137742) q[1];
sx q[1];
rz(-2.9324052) q[1];
sx q[1];
rz(2.7195209) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8764719) q[0];
sx q[0];
rz(-0.19959627) q[0];
sx q[0];
rz(1.100698) q[0];
rz(-1.7447168) q[2];
sx q[2];
rz(-1.2971767) q[2];
sx q[2];
rz(1.5652986) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.6032519) q[1];
sx q[1];
rz(-1.1900239) q[1];
sx q[1];
rz(0.095007665) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.30768577) q[3];
sx q[3];
rz(-0.82594508) q[3];
sx q[3];
rz(0.56203466) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.813039) q[2];
sx q[2];
rz(-2.1812794) q[2];
sx q[2];
rz(-0.73406827) q[2];
rz(-2.0640533) q[3];
sx q[3];
rz(-0.68818337) q[3];
sx q[3];
rz(0.075210007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
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
rz(1.5828736) q[0];
sx q[0];
rz(-2.131077) q[0];
sx q[0];
rz(3.1251113) q[0];
rz(3.0670498) q[1];
sx q[1];
rz(-2.6838979) q[1];
sx q[1];
rz(-2.8864536) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1265565) q[0];
sx q[0];
rz(-1.4269967) q[0];
sx q[0];
rz(2.6565927) q[0];
rz(-1.0696696) q[2];
sx q[2];
rz(-2.5146896) q[2];
sx q[2];
rz(1.3581585) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.8447579) q[1];
sx q[1];
rz(-2.1898309) q[1];
sx q[1];
rz(-2.1059787) q[1];
rz(-pi) q[2];
rz(-1.2665073) q[3];
sx q[3];
rz(-0.36484584) q[3];
sx q[3];
rz(-1.5740652) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.61465803) q[2];
sx q[2];
rz(-2.4444828) q[2];
sx q[2];
rz(-2.441791) q[2];
rz(0.091863306) q[3];
sx q[3];
rz(-1.9242761) q[3];
sx q[3];
rz(2.6868668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5192473) q[0];
sx q[0];
rz(-2.3212101) q[0];
sx q[0];
rz(-2.0297594) q[0];
rz(-0.19019292) q[1];
sx q[1];
rz(-0.88514248) q[1];
sx q[1];
rz(1.9817188) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96944189) q[0];
sx q[0];
rz(-1.5363201) q[0];
sx q[0];
rz(1.2269173) q[0];
rz(1.7669452) q[2];
sx q[2];
rz(-0.64764678) q[2];
sx q[2];
rz(0.37774936) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.1777183) q[1];
sx q[1];
rz(-1.323112) q[1];
sx q[1];
rz(-0.56148333) q[1];
rz(0.044951602) q[3];
sx q[3];
rz(-1.049343) q[3];
sx q[3];
rz(-1.5779881) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.5708892) q[2];
sx q[2];
rz(-0.8231701) q[2];
sx q[2];
rz(2.9111351) q[2];
rz(2.0795836) q[3];
sx q[3];
rz(-1.2758723) q[3];
sx q[3];
rz(1.6331858) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45288169) q[0];
sx q[0];
rz(-1.0588366) q[0];
sx q[0];
rz(1.3558615) q[0];
rz(-0.52976766) q[1];
sx q[1];
rz(-0.7822839) q[1];
sx q[1];
rz(-1.9897602) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.065581948) q[0];
sx q[0];
rz(-1.2154723) q[0];
sx q[0];
rz(-1.6760582) q[0];
x q[1];
rz(-1.0291589) q[2];
sx q[2];
rz(-1.6525606) q[2];
sx q[2];
rz(-2.0643016) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.5024912) q[1];
sx q[1];
rz(-1.7646953) q[1];
sx q[1];
rz(-2.7706233) q[1];
rz(-pi) q[2];
rz(-1.6883259) q[3];
sx q[3];
rz(-1.3580048) q[3];
sx q[3];
rz(-0.26373395) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.4749703) q[2];
sx q[2];
rz(-2.5002561) q[2];
sx q[2];
rz(-2.6744911) q[2];
rz(2.1304255) q[3];
sx q[3];
rz(-0.34365383) q[3];
sx q[3];
rz(-1.3694192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2924627) q[0];
sx q[0];
rz(-0.15501538) q[0];
sx q[0];
rz(0.22468654) q[0];
rz(0.16381964) q[1];
sx q[1];
rz(-0.33912173) q[1];
sx q[1];
rz(-0.64635578) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3098329) q[0];
sx q[0];
rz(-2.0947959) q[0];
sx q[0];
rz(-0.28964596) q[0];
rz(-pi) q[1];
rz(-0.21804131) q[2];
sx q[2];
rz(-1.2550809) q[2];
sx q[2];
rz(1.174508) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.19749459) q[1];
sx q[1];
rz(-2.1131599) q[1];
sx q[1];
rz(-1.6587064) q[1];
rz(2.8180653) q[3];
sx q[3];
rz(-0.40641847) q[3];
sx q[3];
rz(-2.714963) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.15567638) q[2];
sx q[2];
rz(-1.6935657) q[2];
sx q[2];
rz(-2.7316366) q[2];
rz(-2.3112467) q[3];
sx q[3];
rz(-2.1928619) q[3];
sx q[3];
rz(1.4478987) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.603867) q[0];
sx q[0];
rz(-0.69320885) q[0];
sx q[0];
rz(-2.895288) q[0];
rz(-0.82659563) q[1];
sx q[1];
rz(-1.7431424) q[1];
sx q[1];
rz(-2.8776317) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4564698) q[0];
sx q[0];
rz(-1.5635661) q[0];
sx q[0];
rz(-1.475322) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3578431) q[2];
sx q[2];
rz(-1.2211868) q[2];
sx q[2];
rz(2.8658681) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.3166271) q[1];
sx q[1];
rz(-0.74483234) q[1];
sx q[1];
rz(-1.9923575) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.006239) q[3];
sx q[3];
rz(-2.1135984) q[3];
sx q[3];
rz(2.8538728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.10494122) q[2];
sx q[2];
rz(-1.8737326) q[2];
sx q[2];
rz(0.71511739) q[2];
rz(-3.0702843) q[3];
sx q[3];
rz(-1.5569867) q[3];
sx q[3];
rz(0.74688545) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
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
rz(1.0095373) q[0];
sx q[0];
rz(-1.4893463) q[0];
sx q[0];
rz(-2.706053) q[0];
rz(2.9938193) q[1];
sx q[1];
rz(-1.7099893) q[1];
sx q[1];
rz(1.6730283) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3776079) q[0];
sx q[0];
rz(-2.0017636) q[0];
sx q[0];
rz(1.3727643) q[0];
x q[1];
rz(0.69738241) q[2];
sx q[2];
rz(-1.5101523) q[2];
sx q[2];
rz(-1.1866807) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.10994153) q[1];
sx q[1];
rz(-1.1835956) q[1];
sx q[1];
rz(-2.9613858) q[1];
rz(-0.47504039) q[3];
sx q[3];
rz(-0.44774017) q[3];
sx q[3];
rz(2.077075) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.74464166) q[2];
sx q[2];
rz(-1.3775185) q[2];
sx q[2];
rz(-2.2486539) q[2];
rz(1.2785771) q[3];
sx q[3];
rz(-1.9847001) q[3];
sx q[3];
rz(1.4634092) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29499149) q[0];
sx q[0];
rz(-2.8104267) q[0];
sx q[0];
rz(0.60605979) q[0];
rz(-3.1312969) q[1];
sx q[1];
rz(-0.98774397) q[1];
sx q[1];
rz(3.0616679) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.851136) q[0];
sx q[0];
rz(-0.62617597) q[0];
sx q[0];
rz(2.0100223) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.299231) q[2];
sx q[2];
rz(-1.7879221) q[2];
sx q[2];
rz(0.44936839) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.2458107) q[1];
sx q[1];
rz(-1.5485829) q[1];
sx q[1];
rz(0.40708812) q[1];
rz(2.9357443) q[3];
sx q[3];
rz(-1.3449114) q[3];
sx q[3];
rz(-3.1354648) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-3.0004878) q[2];
sx q[2];
rz(-2.2147369) q[2];
sx q[2];
rz(2.1263988) q[2];
rz(-2.5403533) q[3];
sx q[3];
rz(-2.4398118) q[3];
sx q[3];
rz(-1.0465485) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64423185) q[0];
sx q[0];
rz(-1.624122) q[0];
sx q[0];
rz(-2.4540785) q[0];
rz(-2.5524706) q[1];
sx q[1];
rz(-0.67710572) q[1];
sx q[1];
rz(2.6893375) q[1];
rz(1.6670139) q[2];
sx q[2];
rz(-2.0141891) q[2];
sx q[2];
rz(-0.63465848) q[2];
rz(2.3169869) q[3];
sx q[3];
rz(-1.6578703) q[3];
sx q[3];
rz(-0.35005611) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
