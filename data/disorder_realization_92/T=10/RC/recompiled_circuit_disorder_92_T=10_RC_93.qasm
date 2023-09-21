OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.5151514) q[0];
sx q[0];
rz(-0.03349537) q[0];
sx q[0];
rz(1.7749696) q[0];
rz(2.0904436) q[1];
sx q[1];
rz(-1.6519974) q[1];
sx q[1];
rz(2.0096013) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11925919) q[0];
sx q[0];
rz(-2.3046937) q[0];
sx q[0];
rz(1.0975305) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9801324) q[2];
sx q[2];
rz(-1.4070639) q[2];
sx q[2];
rz(1.8979567) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.7992026) q[1];
sx q[1];
rz(-2.4203165) q[1];
sx q[1];
rz(0.93625416) q[1];
rz(-1.2855929) q[3];
sx q[3];
rz(-1.0343026) q[3];
sx q[3];
rz(-0.59629089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.2312317) q[2];
sx q[2];
rz(-1.3277162) q[2];
sx q[2];
rz(1.1532016) q[2];
rz(-0.48405805) q[3];
sx q[3];
rz(-2.3279133) q[3];
sx q[3];
rz(2.6822065) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83950481) q[0];
sx q[0];
rz(-2.8110101) q[0];
sx q[0];
rz(2.6385345) q[0];
rz(-1.5867651) q[1];
sx q[1];
rz(-0.70650548) q[1];
sx q[1];
rz(2.9876626) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.126037) q[0];
sx q[0];
rz(-0.50305788) q[0];
sx q[0];
rz(-1.6004827) q[0];
x q[1];
rz(2.8106861) q[2];
sx q[2];
rz(-0.94793301) q[2];
sx q[2];
rz(1.1922342) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.9736204) q[1];
sx q[1];
rz(-0.64860839) q[1];
sx q[1];
rz(0.82778511) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6318195) q[3];
sx q[3];
rz(-1.5715277) q[3];
sx q[3];
rz(2.7754806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.8872035) q[2];
sx q[2];
rz(-0.35787359) q[2];
sx q[2];
rz(1.3228234) q[2];
rz(-1.6555188) q[3];
sx q[3];
rz(-1.5406939) q[3];
sx q[3];
rz(2.4310908) q[3];
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
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.94674295) q[0];
sx q[0];
rz(-2.0236334) q[0];
sx q[0];
rz(0.90993607) q[0];
rz(0.77727708) q[1];
sx q[1];
rz(-2.3060019) q[1];
sx q[1];
rz(-2.1562703) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18704991) q[0];
sx q[0];
rz(-1.7304725) q[0];
sx q[0];
rz(0.91203441) q[0];
rz(-pi) q[1];
rz(0.91395949) q[2];
sx q[2];
rz(-1.8310391) q[2];
sx q[2];
rz(2.7514806) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.2376155) q[1];
sx q[1];
rz(-1.3217889) q[1];
sx q[1];
rz(1.5140837) q[1];
x q[2];
rz(-0.91089532) q[3];
sx q[3];
rz(-1.0064555) q[3];
sx q[3];
rz(2.6499555) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.2788006) q[2];
sx q[2];
rz(-1.1547487) q[2];
sx q[2];
rz(1.2476236) q[2];
rz(2.6990081) q[3];
sx q[3];
rz(-1.7740039) q[3];
sx q[3];
rz(-0.32143337) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9001532) q[0];
sx q[0];
rz(-1.0192008) q[0];
sx q[0];
rz(2.0181657) q[0];
rz(-2.5627047) q[1];
sx q[1];
rz(-1.6826948) q[1];
sx q[1];
rz(1.3935864) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2508586) q[0];
sx q[0];
rz(-2.3669741) q[0];
sx q[0];
rz(0.61268341) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.34947493) q[2];
sx q[2];
rz(-0.45917837) q[2];
sx q[2];
rz(1.0897204) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.77809282) q[1];
sx q[1];
rz(-2.2908195) q[1];
sx q[1];
rz(-0.6196482) q[1];
rz(1.3136775) q[3];
sx q[3];
rz(-2.3895279) q[3];
sx q[3];
rz(-2.5079692) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.1018155) q[2];
sx q[2];
rz(-1.5856051) q[2];
sx q[2];
rz(-0.0021136443) q[2];
rz(0.56143108) q[3];
sx q[3];
rz(-2.1241472) q[3];
sx q[3];
rz(-0.58825618) q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.026022) q[0];
sx q[0];
rz(-2.0176812) q[0];
sx q[0];
rz(-1.9157238) q[0];
rz(-1.6745802) q[1];
sx q[1];
rz(-1.8672698) q[1];
sx q[1];
rz(-1.7747169) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0923094) q[0];
sx q[0];
rz(-1.7607848) q[0];
sx q[0];
rz(2.7490194) q[0];
x q[1];
rz(-2.3943704) q[2];
sx q[2];
rz(-1.5929654) q[2];
sx q[2];
rz(-1.6037461) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.99202427) q[1];
sx q[1];
rz(-2.4300886) q[1];
sx q[1];
rz(-1.4486538) q[1];
x q[2];
rz(2.6392691) q[3];
sx q[3];
rz(-1.2548496) q[3];
sx q[3];
rz(-1.5497006) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.86429578) q[2];
sx q[2];
rz(-1.0884476) q[2];
sx q[2];
rz(-0.40536353) q[2];
rz(-2.6799485) q[3];
sx q[3];
rz(-2.3159537) q[3];
sx q[3];
rz(1.595114) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6500403) q[0];
sx q[0];
rz(-1.7164282) q[0];
sx q[0];
rz(-0.50338411) q[0];
rz(-0.21884306) q[1];
sx q[1];
rz(-1.8914521) q[1];
sx q[1];
rz(-2.4898081) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9572243) q[0];
sx q[0];
rz(-1.3807218) q[0];
sx q[0];
rz(-2.8977872) q[0];
rz(3.0405365) q[2];
sx q[2];
rz(-2.3172242) q[2];
sx q[2];
rz(1.9161759) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.3411322) q[1];
sx q[1];
rz(-2.9974077) q[1];
sx q[1];
rz(-1.9393117) q[1];
rz(-0.34206335) q[3];
sx q[3];
rz(-0.29933375) q[3];
sx q[3];
rz(0.39810668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.51263222) q[2];
sx q[2];
rz(-1.7444538) q[2];
sx q[2];
rz(2.7499278) q[2];
rz(-2.9351249) q[3];
sx q[3];
rz(-0.73409096) q[3];
sx q[3];
rz(0.68968836) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3702635) q[0];
sx q[0];
rz(-1.6917133) q[0];
sx q[0];
rz(-2.1719334) q[0];
rz(-0.58352739) q[1];
sx q[1];
rz(-2.0136166) q[1];
sx q[1];
rz(3.0113509) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1286344) q[0];
sx q[0];
rz(-0.35612125) q[0];
sx q[0];
rz(-0.14333368) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4728751) q[2];
sx q[2];
rz(-1.0895184) q[2];
sx q[2];
rz(-0.95578335) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.11634593) q[1];
sx q[1];
rz(-1.7882008) q[1];
sx q[1];
rz(2.6245963) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0153766) q[3];
sx q[3];
rz(-2.4820699) q[3];
sx q[3];
rz(-1.5809755) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.51817259) q[2];
sx q[2];
rz(-1.5815846) q[2];
sx q[2];
rz(-2.0557892) q[2];
rz(-3.0893677) q[3];
sx q[3];
rz(-1.4445392) q[3];
sx q[3];
rz(-2.0558555) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6770342) q[0];
sx q[0];
rz(-2.1573986) q[0];
sx q[0];
rz(-1.9751256) q[0];
rz(-0.72558609) q[1];
sx q[1];
rz(-1.2936932) q[1];
sx q[1];
rz(2.3988147) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1449435) q[0];
sx q[0];
rz(-1.7631233) q[0];
sx q[0];
rz(2.3924475) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4648415) q[2];
sx q[2];
rz(-1.9988212) q[2];
sx q[2];
rz(-1.2892436) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.0689773) q[1];
sx q[1];
rz(-1.9779357) q[1];
sx q[1];
rz(-0.8168656) q[1];
rz(-pi) q[2];
rz(0.58917134) q[3];
sx q[3];
rz(-2.1450451) q[3];
sx q[3];
rz(1.9598999) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.62721884) q[2];
sx q[2];
rz(-1.5292995) q[2];
sx q[2];
rz(0.26088866) q[2];
rz(-1.1076814) q[3];
sx q[3];
rz(-1.2872144) q[3];
sx q[3];
rz(1.0816983) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35659197) q[0];
sx q[0];
rz(-2.3605425) q[0];
sx q[0];
rz(2.2055431) q[0];
rz(-3.127457) q[1];
sx q[1];
rz(-1.8363258) q[1];
sx q[1];
rz(-2.4900808) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0148841) q[0];
sx q[0];
rz(-1.5471022) q[0];
sx q[0];
rz(3.1372877) q[0];
rz(0.94349761) q[2];
sx q[2];
rz(-2.128696) q[2];
sx q[2];
rz(-0.055659143) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.1956605) q[1];
sx q[1];
rz(-1.4504499) q[1];
sx q[1];
rz(1.2502928) q[1];
rz(-pi) q[2];
rz(-1.2048079) q[3];
sx q[3];
rz(-1.9848739) q[3];
sx q[3];
rz(-0.28596349) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.8462048) q[2];
sx q[2];
rz(-0.69869763) q[2];
sx q[2];
rz(-0.52337581) q[2];
rz(2.8335617) q[3];
sx q[3];
rz(-0.57294661) q[3];
sx q[3];
rz(1.3380922) q[3];
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
rz(-1.4073407) q[0];
sx q[0];
rz(-2.9946406) q[0];
sx q[0];
rz(-1.403632) q[0];
rz(-0.47406468) q[1];
sx q[1];
rz(-1.4539376) q[1];
sx q[1];
rz(1.7636991) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5705469) q[0];
sx q[0];
rz(-2.018398) q[0];
sx q[0];
rz(-0.96121995) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4268742) q[2];
sx q[2];
rz(-0.93313365) q[2];
sx q[2];
rz(-1.4710466) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.5046895) q[1];
sx q[1];
rz(-2.6733589) q[1];
sx q[1];
rz(0.56028985) q[1];
x q[2];
rz(0.5573427) q[3];
sx q[3];
rz(-1.4321657) q[3];
sx q[3];
rz(-2.0687452) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.3906117) q[2];
sx q[2];
rz(-1.5079974) q[2];
sx q[2];
rz(-1.5314468) q[2];
rz(-1.6837998) q[3];
sx q[3];
rz(-2.0861574) q[3];
sx q[3];
rz(-2.2916601) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67125852) q[0];
sx q[0];
rz(-0.27161921) q[0];
sx q[0];
rz(-2.0441396) q[0];
rz(0.63256565) q[1];
sx q[1];
rz(-1.0610776) q[1];
sx q[1];
rz(-2.9656596) q[1];
rz(1.8995646) q[2];
sx q[2];
rz(-2.2148) q[2];
sx q[2];
rz(2.7684341) q[2];
rz(-2.408151) q[3];
sx q[3];
rz(-2.1248795) q[3];
sx q[3];
rz(-2.636573) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
