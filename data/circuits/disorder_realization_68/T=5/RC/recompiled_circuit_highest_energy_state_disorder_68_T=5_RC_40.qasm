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
rz(-2.3013378) q[0];
sx q[0];
rz(-1.037896) q[0];
sx q[0];
rz(2.2965746) q[0];
rz(-2.4294699) q[1];
sx q[1];
rz(-1.0016088) q[1];
sx q[1];
rz(-1.4955624) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9121542) q[0];
sx q[0];
rz(-0.99675465) q[0];
sx q[0];
rz(-2.6340934) q[0];
rz(-pi) q[1];
rz(0.38365429) q[2];
sx q[2];
rz(-1.7100514) q[2];
sx q[2];
rz(2.4568707) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.6604267) q[1];
sx q[1];
rz(-1.1200953) q[1];
sx q[1];
rz(-0.53498241) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5450685) q[3];
sx q[3];
rz(-0.46677315) q[3];
sx q[3];
rz(2.1201565) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.738395) q[2];
sx q[2];
rz(-1.1824111) q[2];
sx q[2];
rz(2.5851868) q[2];
rz(2.6465936) q[3];
sx q[3];
rz(-0.35111108) q[3];
sx q[3];
rz(0.88465148) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2719088) q[0];
sx q[0];
rz(-0.10232919) q[0];
sx q[0];
rz(2.5129357) q[0];
rz(-2.2643845) q[1];
sx q[1];
rz(-2.6429206) q[1];
sx q[1];
rz(0.25310755) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0639187) q[0];
sx q[0];
rz(-1.5787933) q[0];
sx q[0];
rz(-1.5863938) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.44127522) q[2];
sx q[2];
rz(-1.8417449) q[2];
sx q[2];
rz(0.51368062) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7444708) q[1];
sx q[1];
rz(-1.4543449) q[1];
sx q[1];
rz(2.9411903) q[1];
rz(-pi) q[2];
x q[2];
rz(0.25365003) q[3];
sx q[3];
rz(-2.2861655) q[3];
sx q[3];
rz(-0.77798346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.46193281) q[2];
sx q[2];
rz(-1.9340645) q[2];
sx q[2];
rz(1.0665464) q[2];
rz(1.3072394) q[3];
sx q[3];
rz(-2.6756838) q[3];
sx q[3];
rz(0.90827847) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9268554) q[0];
sx q[0];
rz(-2.1934788) q[0];
sx q[0];
rz(-0.98534775) q[0];
rz(2.7477879) q[1];
sx q[1];
rz(-2.1486798) q[1];
sx q[1];
rz(2.6148112) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5827127) q[0];
sx q[0];
rz(-1.0436907) q[0];
sx q[0];
rz(0.096313535) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5622732) q[2];
sx q[2];
rz(-2.0181106) q[2];
sx q[2];
rz(-1.2944702) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.56934988) q[1];
sx q[1];
rz(-1.7399638) q[1];
sx q[1];
rz(1.1430686) q[1];
rz(-pi) q[2];
rz(2.18264) q[3];
sx q[3];
rz(-0.69330319) q[3];
sx q[3];
rz(0.43402754) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.3452611) q[2];
sx q[2];
rz(-0.76097208) q[2];
sx q[2];
rz(-2.995028) q[2];
rz(-0.86756271) q[3];
sx q[3];
rz(-2.2680794) q[3];
sx q[3];
rz(2.383702) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4124311) q[0];
sx q[0];
rz(-0.49598345) q[0];
sx q[0];
rz(2.5878986) q[0];
rz(-1.6665392) q[1];
sx q[1];
rz(-0.29860425) q[1];
sx q[1];
rz(0.24436229) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9825443) q[0];
sx q[0];
rz(-1.9031018) q[0];
sx q[0];
rz(-0.37489069) q[0];
rz(-pi) q[1];
rz(2.8344523) q[2];
sx q[2];
rz(-2.2835586) q[2];
sx q[2];
rz(-2.3736726) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.7556954) q[1];
sx q[1];
rz(-1.9287957) q[1];
sx q[1];
rz(2.6935607) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5457754) q[3];
sx q[3];
rz(-0.85880781) q[3];
sx q[3];
rz(2.156949) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.060140572) q[2];
sx q[2];
rz(-1.1271366) q[2];
sx q[2];
rz(-0.78090182) q[2];
rz(2.7447356) q[3];
sx q[3];
rz(-3.1271827) q[3];
sx q[3];
rz(-1.074033) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2395372) q[0];
sx q[0];
rz(-0.1665512) q[0];
sx q[0];
rz(0.2567513) q[0];
rz(2.9638839) q[1];
sx q[1];
rz(-0.54568988) q[1];
sx q[1];
rz(0.94388747) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2985757) q[0];
sx q[0];
rz(-2.0110879) q[0];
sx q[0];
rz(0.40233516) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2805528) q[2];
sx q[2];
rz(-1.3849045) q[2];
sx q[2];
rz(-3.0772532) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.1348995) q[1];
sx q[1];
rz(-2.1104782) q[1];
sx q[1];
rz(1.5758908) q[1];
x q[2];
rz(2.2511399) q[3];
sx q[3];
rz(-0.6750921) q[3];
sx q[3];
rz(-2.4321041) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.2029767) q[2];
sx q[2];
rz(-1.0304893) q[2];
sx q[2];
rz(-0.20982783) q[2];
rz(-3.0948011) q[3];
sx q[3];
rz(-1.4593461) q[3];
sx q[3];
rz(-0.34745026) q[3];
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
rz(-0.059134722) q[0];
sx q[0];
rz(-2.3248398) q[0];
sx q[0];
rz(-1.2030075) q[0];
rz(1.5164392) q[1];
sx q[1];
rz(-2.7639183) q[1];
sx q[1];
rz(0.032940544) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0982978) q[0];
sx q[0];
rz(-1.2970719) q[0];
sx q[0];
rz(2.0541463) q[0];
rz(-pi) q[1];
rz(0.85196544) q[2];
sx q[2];
rz(-0.72418606) q[2];
sx q[2];
rz(2.9018096) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.6194544) q[1];
sx q[1];
rz(-0.65548766) q[1];
sx q[1];
rz(-0.54490276) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.88992702) q[3];
sx q[3];
rz(-1.0491199) q[3];
sx q[3];
rz(1.8311178) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.5861627) q[2];
sx q[2];
rz(-2.1383492) q[2];
sx q[2];
rz(0.45247751) q[2];
rz(0.14719506) q[3];
sx q[3];
rz(-0.23284027) q[3];
sx q[3];
rz(1.8424621) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.074987) q[0];
sx q[0];
rz(-0.43656483) q[0];
sx q[0];
rz(2.7845352) q[0];
rz(-1.0128516) q[1];
sx q[1];
rz(-1.9722936) q[1];
sx q[1];
rz(3.1049407) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9215611) q[0];
sx q[0];
rz(-2.4274094) q[0];
sx q[0];
rz(-2.4490687) q[0];
rz(1.4674856) q[2];
sx q[2];
rz(-1.1432308) q[2];
sx q[2];
rz(2.4201938) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.8690435) q[1];
sx q[1];
rz(-0.70768917) q[1];
sx q[1];
rz(-1.7131192) q[1];
rz(-0.030825214) q[3];
sx q[3];
rz(-1.5753058) q[3];
sx q[3];
rz(-2.4140178) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.4770294) q[2];
sx q[2];
rz(-0.9534812) q[2];
sx q[2];
rz(-3.0518517) q[2];
rz(-1.2612032) q[3];
sx q[3];
rz(-1.829105) q[3];
sx q[3];
rz(1.7173654) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4795714) q[0];
sx q[0];
rz(-2.5630072) q[0];
sx q[0];
rz(1.660996) q[0];
rz(1.6914233) q[1];
sx q[1];
rz(-0.65933508) q[1];
sx q[1];
rz(-0.75993842) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63045365) q[0];
sx q[0];
rz(-1.5811111) q[0];
sx q[0];
rz(1.4546088) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.12822656) q[2];
sx q[2];
rz(-2.5094446) q[2];
sx q[2];
rz(1.8116443) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.6531585) q[1];
sx q[1];
rz(-1.9699082) q[1];
sx q[1];
rz(2.0293971) q[1];
rz(-pi) q[2];
rz(1.8283056) q[3];
sx q[3];
rz(-2.0638568) q[3];
sx q[3];
rz(2.8185533) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.30457589) q[2];
sx q[2];
rz(-1.872007) q[2];
sx q[2];
rz(0.059583511) q[2];
rz(0.42851055) q[3];
sx q[3];
rz(-0.84463745) q[3];
sx q[3];
rz(2.5394411) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8932327) q[0];
sx q[0];
rz(-0.28286523) q[0];
sx q[0];
rz(-2.7981113) q[0];
rz(0.19613014) q[1];
sx q[1];
rz(-0.88534147) q[1];
sx q[1];
rz(0.79998618) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72092743) q[0];
sx q[0];
rz(-1.979083) q[0];
sx q[0];
rz(-0.83939206) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3197199) q[2];
sx q[2];
rz(-1.0227239) q[2];
sx q[2];
rz(2.1679116) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.6226253) q[1];
sx q[1];
rz(-2.1299348) q[1];
sx q[1];
rz(-1.6230727) q[1];
x q[2];
rz(-3.0723582) q[3];
sx q[3];
rz(-0.88919176) q[3];
sx q[3];
rz(0.72661663) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.7949152) q[2];
sx q[2];
rz(-0.76378834) q[2];
sx q[2];
rz(-2.9128892) q[2];
rz(-1.4165357) q[3];
sx q[3];
rz(-1.4226457) q[3];
sx q[3];
rz(-2.4712839) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5423841) q[0];
sx q[0];
rz(-0.51775652) q[0];
sx q[0];
rz(1.9589348) q[0];
rz(-2.5987437) q[1];
sx q[1];
rz(-0.12195568) q[1];
sx q[1];
rz(0.45546946) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1802954) q[0];
sx q[0];
rz(-2.0019128) q[0];
sx q[0];
rz(-2.6102553) q[0];
rz(-1.9827911) q[2];
sx q[2];
rz(-0.48589009) q[2];
sx q[2];
rz(2.0249572) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.7389981) q[1];
sx q[1];
rz(-1.7070424) q[1];
sx q[1];
rz(2.6722277) q[1];
rz(-pi) q[2];
rz(0.72584589) q[3];
sx q[3];
rz(-0.56768394) q[3];
sx q[3];
rz(0.7782225) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.9246284) q[2];
sx q[2];
rz(-2.2483726) q[2];
sx q[2];
rz(-2.7157937) q[2];
rz(0.18481542) q[3];
sx q[3];
rz(-0.63822377) q[3];
sx q[3];
rz(3.1137915) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4375147) q[0];
sx q[0];
rz(-1.5298433) q[0];
sx q[0];
rz(-0.035506305) q[0];
rz(-0.47954814) q[1];
sx q[1];
rz(-1.5621114) q[1];
sx q[1];
rz(1.5893804) q[1];
rz(0.32301474) q[2];
sx q[2];
rz(-2.5217006) q[2];
sx q[2];
rz(2.9519242) q[2];
rz(-1.7749589) q[3];
sx q[3];
rz(-0.62075172) q[3];
sx q[3];
rz(-0.78637607) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
