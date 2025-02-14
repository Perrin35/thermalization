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
rz(-0.79987502) q[0];
sx q[0];
rz(2.3246111) q[0];
sx q[0];
rz(9.8819879) q[0];
rz(0.41181052) q[1];
sx q[1];
rz(-1.5195941) q[1];
sx q[1];
rz(-2.5817459) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7871817) q[0];
sx q[0];
rz(-2.5776064) q[0];
sx q[0];
rz(0.63040479) q[0];
rz(2.0617606) q[2];
sx q[2];
rz(-0.55779558) q[2];
sx q[2];
rz(1.3809134) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.76418855) q[1];
sx q[1];
rz(-1.848583) q[1];
sx q[1];
rz(2.6242816) q[1];
rz(-pi) q[2];
rz(-2.9212037) q[3];
sx q[3];
rz(-1.427009) q[3];
sx q[3];
rz(1.4291562) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.7655699) q[2];
sx q[2];
rz(-0.96161181) q[2];
sx q[2];
rz(1.8454856) q[2];
rz(0.62093312) q[3];
sx q[3];
rz(-0.66578484) q[3];
sx q[3];
rz(-0.92500979) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6821482) q[0];
sx q[0];
rz(-2.5542673) q[0];
sx q[0];
rz(-2.4272954) q[0];
rz(-0.722305) q[1];
sx q[1];
rz(-1.0853094) q[1];
sx q[1];
rz(1.8169656) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9324397) q[0];
sx q[0];
rz(-2.5896833) q[0];
sx q[0];
rz(0.96262424) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5157976) q[2];
sx q[2];
rz(-2.3028513) q[2];
sx q[2];
rz(-2.3961803) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.2683619) q[1];
sx q[1];
rz(-1.2030081) q[1];
sx q[1];
rz(1.7413543) q[1];
x q[2];
rz(1.7198661) q[3];
sx q[3];
rz(-2.4507634) q[3];
sx q[3];
rz(-0.2836424) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.033919949) q[2];
sx q[2];
rz(-1.718797) q[2];
sx q[2];
rz(-0.99204341) q[2];
rz(-0.77573675) q[3];
sx q[3];
rz(-2.3596767) q[3];
sx q[3];
rz(-1.0824664) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42763212) q[0];
sx q[0];
rz(-1.3466703) q[0];
sx q[0];
rz(-2.2295075) q[0];
rz(-0.38463792) q[1];
sx q[1];
rz(-0.96133989) q[1];
sx q[1];
rz(-2.6461163) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3289483) q[0];
sx q[0];
rz(-1.5248796) q[0];
sx q[0];
rz(-1.4981235) q[0];
rz(-pi) q[1];
rz(-0.20553164) q[2];
sx q[2];
rz(-2.3816817) q[2];
sx q[2];
rz(0.854137) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.8532995) q[1];
sx q[1];
rz(-1.5202731) q[1];
sx q[1];
rz(-2.4171197) q[1];
rz(-pi) q[2];
x q[2];
rz(1.552565) q[3];
sx q[3];
rz(-2.5196001) q[3];
sx q[3];
rz(2.7014159) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.89874011) q[2];
sx q[2];
rz(-2.0169368) q[2];
sx q[2];
rz(2.7070572) q[2];
rz(3.0430326) q[3];
sx q[3];
rz(-1.7193272) q[3];
sx q[3];
rz(0.77384531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
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
rz(-0.74715215) q[0];
sx q[0];
rz(-0.23500615) q[0];
sx q[0];
rz(2.8727942) q[0];
rz(-1.0768184) q[1];
sx q[1];
rz(-1.7348758) q[1];
sx q[1];
rz(1.1928308) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6859602) q[0];
sx q[0];
rz(-1.3979027) q[0];
sx q[0];
rz(0.25211035) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4401764) q[2];
sx q[2];
rz(-1.7057749) q[2];
sx q[2];
rz(-0.93739742) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.12700203) q[1];
sx q[1];
rz(-1.4140097) q[1];
sx q[1];
rz(0.21611045) q[1];
rz(-pi) q[2];
rz(2.4402839) q[3];
sx q[3];
rz(-1.7265254) q[3];
sx q[3];
rz(2.6913672) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.8863135) q[2];
sx q[2];
rz(-2.2463319) q[2];
sx q[2];
rz(0.2992343) q[2];
rz(2.8507161) q[3];
sx q[3];
rz(-1.2521005) q[3];
sx q[3];
rz(0.65557426) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3000779) q[0];
sx q[0];
rz(-2.1617007) q[0];
sx q[0];
rz(3.063524) q[0];
rz(2.1878751) q[1];
sx q[1];
rz(-2.6225312) q[1];
sx q[1];
rz(1.5740707) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.038356) q[0];
sx q[0];
rz(-0.7013196) q[0];
sx q[0];
rz(-1.288802) q[0];
rz(-pi) q[1];
rz(1.2260796) q[2];
sx q[2];
rz(-0.90391747) q[2];
sx q[2];
rz(2.4645811) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5213274) q[1];
sx q[1];
rz(-0.6129188) q[1];
sx q[1];
rz(1.6271724) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5007581) q[3];
sx q[3];
rz(-0.87890714) q[3];
sx q[3];
rz(2.6393011) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.8629525) q[2];
sx q[2];
rz(-0.43206698) q[2];
sx q[2];
rz(-2.8734109) q[2];
rz(-2.6523318) q[3];
sx q[3];
rz(-2.0475976) q[3];
sx q[3];
rz(-1.6930273) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(-2.0972524) q[0];
sx q[0];
rz(-3.1179929) q[0];
sx q[0];
rz(-2.6771255) q[0];
rz(2.6181472) q[1];
sx q[1];
rz(-2.372066) q[1];
sx q[1];
rz(-0.73062599) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0574065) q[0];
sx q[0];
rz(-2.3825186) q[0];
sx q[0];
rz(1.7763863) q[0];
rz(-pi) q[1];
rz(-2.7774657) q[2];
sx q[2];
rz(-2.4362323) q[2];
sx q[2];
rz(-3.0400227) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.1878512) q[1];
sx q[1];
rz(-1.9914728) q[1];
sx q[1];
rz(-0.71297808) q[1];
x q[2];
rz(-1.4387801) q[3];
sx q[3];
rz(-1.8204597) q[3];
sx q[3];
rz(0.77900797) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.0475433) q[2];
sx q[2];
rz(-2.0619679) q[2];
sx q[2];
rz(3.003982) q[2];
rz(-1.1329457) q[3];
sx q[3];
rz(-1.1762985) q[3];
sx q[3];
rz(1.2631811) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15261821) q[0];
sx q[0];
rz(-1.512383) q[0];
sx q[0];
rz(-0.033893943) q[0];
rz(-2.7821817) q[1];
sx q[1];
rz(-1.5243328) q[1];
sx q[1];
rz(0.83438897) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8501668) q[0];
sx q[0];
rz(-1.8803673) q[0];
sx q[0];
rz(-0.39637027) q[0];
rz(-pi) q[1];
rz(-2.5691146) q[2];
sx q[2];
rz(-2.3668063) q[2];
sx q[2];
rz(3.0511193) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5339116) q[1];
sx q[1];
rz(-1.6661196) q[1];
sx q[1];
rz(-1.3168112) q[1];
rz(-2.4700048) q[3];
sx q[3];
rz(-1.62) q[3];
sx q[3];
rz(-2.4528596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.413201) q[2];
sx q[2];
rz(-0.30343702) q[2];
sx q[2];
rz(-1.0580074) q[2];
rz(-0.47438619) q[3];
sx q[3];
rz(-2.3460903) q[3];
sx q[3];
rz(2.0856196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0673085) q[0];
sx q[0];
rz(-1.625076) q[0];
sx q[0];
rz(-1.7838595) q[0];
rz(2.7108497) q[1];
sx q[1];
rz(-1.4746702) q[1];
sx q[1];
rz(-2.6745083) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3245832) q[0];
sx q[0];
rz(-1.5276434) q[0];
sx q[0];
rz(0.0014716455) q[0];
x q[1];
rz(1.378304) q[2];
sx q[2];
rz(-0.87710947) q[2];
sx q[2];
rz(0.81499962) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.7469646) q[1];
sx q[1];
rz(-2.10022) q[1];
sx q[1];
rz(-1.5922597) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0881431) q[3];
sx q[3];
rz(-0.9336144) q[3];
sx q[3];
rz(0.11252014) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0908541) q[2];
sx q[2];
rz(-2.723912) q[2];
sx q[2];
rz(2.0609071) q[2];
rz(2.4596227) q[3];
sx q[3];
rz(-2.3365648) q[3];
sx q[3];
rz(0.69847703) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4526116) q[0];
sx q[0];
rz(-0.51849759) q[0];
sx q[0];
rz(2.3338351) q[0];
rz(-2.141433) q[1];
sx q[1];
rz(-0.88614416) q[1];
sx q[1];
rz(-2.3449786) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33045443) q[0];
sx q[0];
rz(-1.66798) q[0];
sx q[0];
rz(-3.099346) q[0];
rz(0.62532897) q[2];
sx q[2];
rz(-1.3742067) q[2];
sx q[2];
rz(-2.1670053) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.63846) q[1];
sx q[1];
rz(-1.1640932) q[1];
sx q[1];
rz(-0.15127123) q[1];
rz(-1.2563989) q[3];
sx q[3];
rz(-1.8125497) q[3];
sx q[3];
rz(0.35255656) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.5920068) q[2];
sx q[2];
rz(-1.0909811) q[2];
sx q[2];
rz(-0.31527147) q[2];
rz(-1.573805) q[3];
sx q[3];
rz(-1.9939634) q[3];
sx q[3];
rz(-2.018759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0386117) q[0];
sx q[0];
rz(-2.4760315) q[0];
sx q[0];
rz(0.82157201) q[0];
rz(-2.6577677) q[1];
sx q[1];
rz(-1.7211823) q[1];
sx q[1];
rz(-0.22463591) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.076775) q[0];
sx q[0];
rz(-0.358264) q[0];
sx q[0];
rz(-1.129877) q[0];
rz(-pi) q[1];
x q[1];
rz(0.63172747) q[2];
sx q[2];
rz(-1.016482) q[2];
sx q[2];
rz(0.71507031) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.37525081) q[1];
sx q[1];
rz(-0.90000376) q[1];
sx q[1];
rz(0.2672345) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7849237) q[3];
sx q[3];
rz(-1.4487293) q[3];
sx q[3];
rz(-1.7904953) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.75767526) q[2];
sx q[2];
rz(-3.0249247) q[2];
sx q[2];
rz(-0.095001027) q[2];
rz(1.4875686) q[3];
sx q[3];
rz(-2.5520958) q[3];
sx q[3];
rz(2.3768363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8944396) q[0];
sx q[0];
rz(-2.7338487) q[0];
sx q[0];
rz(-0.26640531) q[0];
rz(2.5447625) q[1];
sx q[1];
rz(-1.4529556) q[1];
sx q[1];
rz(-1.6642889) q[1];
rz(2.785037) q[2];
sx q[2];
rz(-2.2168163) q[2];
sx q[2];
rz(1.4067895) q[2];
rz(1.5726907) q[3];
sx q[3];
rz(-1.9463149) q[3];
sx q[3];
rz(-2.7494242) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
