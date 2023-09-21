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
rz(-1.3666231) q[0];
rz(2.0904436) q[1];
sx q[1];
rz(-1.6519974) q[1];
sx q[1];
rz(-1.1319914) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6080947) q[0];
sx q[0];
rz(-0.84871263) q[0];
sx q[0];
rz(-0.46790926) q[0];
rz(-2.3425383) q[2];
sx q[2];
rz(-0.22944268) q[2];
sx q[2];
rz(0.45861751) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.0238029) q[1];
sx q[1];
rz(-2.1315247) q[1];
sx q[1];
rz(0.48052131) q[1];
x q[2];
rz(-2.6996783) q[3];
sx q[3];
rz(-0.60097296) q[3];
sx q[3];
rz(1.1170944) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.91036096) q[2];
sx q[2];
rz(-1.3277162) q[2];
sx q[2];
rz(-1.1532016) q[2];
rz(0.48405805) q[3];
sx q[3];
rz(-0.81367937) q[3];
sx q[3];
rz(2.6822065) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83950481) q[0];
sx q[0];
rz(-0.33058259) q[0];
sx q[0];
rz(-0.50305811) q[0];
rz(1.5867651) q[1];
sx q[1];
rz(-2.4350872) q[1];
sx q[1];
rz(2.9876626) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0155556) q[0];
sx q[0];
rz(-2.6385348) q[0];
sx q[0];
rz(1.6004827) q[0];
x q[1];
rz(0.33090654) q[2];
sx q[2];
rz(-2.1936596) q[2];
sx q[2];
rz(1.1922342) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.0347552) q[1];
sx q[1];
rz(-1.1498067) q[1];
sx q[1];
rz(1.0616598) q[1];
rz(-pi) q[2];
rz(3.1400938) q[3];
sx q[3];
rz(-2.631819) q[3];
sx q[3];
rz(-1.203376) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.2543891) q[2];
sx q[2];
rz(-2.7837191) q[2];
sx q[2];
rz(1.8187693) q[2];
rz(1.6555188) q[3];
sx q[3];
rz(-1.5406939) q[3];
sx q[3];
rz(-2.4310908) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94674295) q[0];
sx q[0];
rz(-2.0236334) q[0];
sx q[0];
rz(-0.90993607) q[0];
rz(-2.3643156) q[1];
sx q[1];
rz(-0.83559075) q[1];
sx q[1];
rz(2.1562703) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9604208) q[0];
sx q[0];
rz(-0.67502484) q[0];
sx q[0];
rz(1.8280562) q[0];
rz(-pi) q[1];
rz(1.9820205) q[2];
sx q[2];
rz(-0.69934884) q[2];
sx q[2];
rz(-1.6388091) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.2376155) q[1];
sx q[1];
rz(-1.3217889) q[1];
sx q[1];
rz(-1.5140837) q[1];
rz(2.3722234) q[3];
sx q[3];
rz(-2.3017075) q[3];
sx q[3];
rz(-1.6826671) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.862792) q[2];
sx q[2];
rz(-1.9868439) q[2];
sx q[2];
rz(1.8939691) q[2];
rz(0.4425846) q[3];
sx q[3];
rz(-1.3675888) q[3];
sx q[3];
rz(2.8201593) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9001532) q[0];
sx q[0];
rz(-2.1223919) q[0];
sx q[0];
rz(1.1234269) q[0];
rz(-2.5627047) q[1];
sx q[1];
rz(-1.4588979) q[1];
sx q[1];
rz(1.7480063) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9959895) q[0];
sx q[0];
rz(-1.1568501) q[0];
sx q[0];
rz(2.4664509) q[0];
x q[1];
rz(-0.43487866) q[2];
sx q[2];
rz(-1.41845) q[2];
sx q[2];
rz(-0.16532126) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.046603831) q[1];
sx q[1];
rz(-0.91218439) q[1];
sx q[1];
rz(-0.98595001) q[1];
rz(0.23354236) q[3];
sx q[3];
rz(-0.84905784) q[3];
sx q[3];
rz(2.8535709) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.1018155) q[2];
sx q[2];
rz(-1.5559876) q[2];
sx q[2];
rz(-0.0021136443) q[2];
rz(2.5801616) q[3];
sx q[3];
rz(-1.0174454) q[3];
sx q[3];
rz(-0.58825618) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11557065) q[0];
sx q[0];
rz(-2.0176812) q[0];
sx q[0];
rz(-1.9157238) q[0];
rz(-1.6745802) q[1];
sx q[1];
rz(-1.2743228) q[1];
sx q[1];
rz(1.7747169) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0492832) q[0];
sx q[0];
rz(-1.7607848) q[0];
sx q[0];
rz(-2.7490194) q[0];
x q[1];
rz(0.032614313) q[2];
sx q[2];
rz(-2.3941052) q[2];
sx q[2];
rz(3.0847197) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.83134507) q[1];
sx q[1];
rz(-2.275895) q[1];
sx q[1];
rz(-3.0369333) q[1];
rz(-pi) q[2];
rz(2.5451238) q[3];
sx q[3];
rz(-0.58613741) q[3];
sx q[3];
rz(-0.49367192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.86429578) q[2];
sx q[2];
rz(-1.0884476) q[2];
sx q[2];
rz(-2.7362291) q[2];
rz(-2.6799485) q[3];
sx q[3];
rz(-2.3159537) q[3];
sx q[3];
rz(1.595114) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49155238) q[0];
sx q[0];
rz(-1.4251645) q[0];
sx q[0];
rz(-2.6382085) q[0];
rz(0.21884306) q[1];
sx q[1];
rz(-1.8914521) q[1];
sx q[1];
rz(-0.65178451) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1056021) q[0];
sx q[0];
rz(-2.833617) q[0];
sx q[0];
rz(-0.67291798) q[0];
rz(3.0405365) q[2];
sx q[2];
rz(-0.82436845) q[2];
sx q[2];
rz(-1.9161759) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.5946878) q[1];
sx q[1];
rz(-1.5190131) q[1];
sx q[1];
rz(1.7054218) q[1];
x q[2];
rz(0.34206335) q[3];
sx q[3];
rz(-0.29933375) q[3];
sx q[3];
rz(-0.39810668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.6289604) q[2];
sx q[2];
rz(-1.3971389) q[2];
sx q[2];
rz(-0.39166489) q[2];
rz(-0.20646778) q[3];
sx q[3];
rz(-0.73409096) q[3];
sx q[3];
rz(2.4519043) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3702635) q[0];
sx q[0];
rz(-1.6917133) q[0];
sx q[0];
rz(-0.96965924) q[0];
rz(-2.5580653) q[1];
sx q[1];
rz(-2.0136166) q[1];
sx q[1];
rz(-3.0113509) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1657432) q[0];
sx q[0];
rz(-1.2184869) q[0];
sx q[0];
rz(-1.5177112) q[0];
x q[1];
rz(0.98353705) q[2];
sx q[2];
rz(-2.1526255) q[2];
sx q[2];
rz(2.8772417) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.8091734) q[1];
sx q[1];
rz(-2.0744588) q[1];
sx q[1];
rz(1.8196351) q[1];
rz(2.7534915) q[3];
sx q[3];
rz(-1.023205) q[3];
sx q[3];
rz(-2.2263118) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.51817259) q[2];
sx q[2];
rz(-1.5815846) q[2];
sx q[2];
rz(-1.0858034) q[2];
rz(3.0893677) q[3];
sx q[3];
rz(-1.4445392) q[3];
sx q[3];
rz(2.0558555) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4645585) q[0];
sx q[0];
rz(-2.1573986) q[0];
sx q[0];
rz(1.1664671) q[0];
rz(-2.4160066) q[1];
sx q[1];
rz(-1.8478994) q[1];
sx q[1];
rz(-0.74277791) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.60177875) q[0];
sx q[0];
rz(-0.83866461) q[0];
sx q[0];
rz(1.3108805) q[0];
rz(-2.5119945) q[2];
sx q[2];
rz(-0.78231914) q[2];
sx q[2];
rz(-2.9462189) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.0689773) q[1];
sx q[1];
rz(-1.9779357) q[1];
sx q[1];
rz(2.3247271) q[1];
x q[2];
rz(-0.58917134) q[3];
sx q[3];
rz(-0.99654752) q[3];
sx q[3];
rz(1.9598999) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.5143738) q[2];
sx q[2];
rz(-1.5292995) q[2];
sx q[2];
rz(0.26088866) q[2];
rz(2.0339113) q[3];
sx q[3];
rz(-1.8543782) q[3];
sx q[3];
rz(-1.0816983) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35659197) q[0];
sx q[0];
rz(-0.78105015) q[0];
sx q[0];
rz(2.2055431) q[0];
rz(-0.014135663) q[1];
sx q[1];
rz(-1.8363258) q[1];
sx q[1];
rz(-0.65151185) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0148841) q[0];
sx q[0];
rz(-1.5471022) q[0];
sx q[0];
rz(-3.1372877) q[0];
rz(-0.94349761) q[2];
sx q[2];
rz(-1.0128967) q[2];
sx q[2];
rz(3.0859335) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.028101746) q[1];
sx q[1];
rz(-2.7999702) q[1];
sx q[1];
rz(-1.9373059) q[1];
rz(-pi) q[2];
rz(-2.7016919) q[3];
sx q[3];
rz(-1.2370046) q[3];
sx q[3];
rz(-1.7037638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.8462048) q[2];
sx q[2];
rz(-2.442895) q[2];
sx q[2];
rz(-2.6182168) q[2];
rz(-0.30803099) q[3];
sx q[3];
rz(-0.57294661) q[3];
sx q[3];
rz(1.3380922) q[3];
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
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.734252) q[0];
sx q[0];
rz(-2.9946406) q[0];
sx q[0];
rz(-1.7379606) q[0];
rz(-0.47406468) q[1];
sx q[1];
rz(-1.4539376) q[1];
sx q[1];
rz(-1.3778936) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70625593) q[0];
sx q[0];
rz(-1.0284817) q[0];
sx q[0];
rz(-0.52973912) q[0];
rz(2.2950298) q[2];
sx q[2];
rz(-0.91869527) q[2];
sx q[2];
rz(-0.50154274) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-3.1174406) q[1];
sx q[1];
rz(-1.9630868) q[1];
sx q[1];
rz(1.8333608) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4078478) q[3];
sx q[3];
rz(-2.1221707) q[3];
sx q[3];
rz(0.58386246) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.3906117) q[2];
sx q[2];
rz(-1.6335952) q[2];
sx q[2];
rz(-1.5314468) q[2];
rz(-1.4577929) q[3];
sx q[3];
rz(-1.0554353) q[3];
sx q[3];
rz(0.84993258) q[3];
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
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67125852) q[0];
sx q[0];
rz(-0.27161921) q[0];
sx q[0];
rz(-2.0441396) q[0];
rz(2.509027) q[1];
sx q[1];
rz(-2.0805151) q[1];
sx q[1];
rz(0.17593304) q[1];
rz(2.4709672) q[2];
sx q[2];
rz(-1.8319596) q[2];
sx q[2];
rz(-1.7419227) q[2];
rz(-2.395527) q[3];
sx q[3];
rz(-2.2545771) q[3];
sx q[3];
rz(1.547326) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
