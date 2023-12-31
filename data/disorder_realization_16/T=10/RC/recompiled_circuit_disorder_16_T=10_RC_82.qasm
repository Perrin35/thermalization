OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.47473946) q[0];
sx q[0];
rz(-0.82959509) q[0];
sx q[0];
rz(0.15396804) q[0];
rz(-2.3078168) q[1];
sx q[1];
rz(-0.99234617) q[1];
sx q[1];
rz(0.33831236) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6217123) q[0];
sx q[0];
rz(-1.3268688) q[0];
sx q[0];
rz(1.8983311) q[0];
rz(-pi) q[1];
rz(-2.0182274) q[2];
sx q[2];
rz(-0.72927232) q[2];
sx q[2];
rz(0.73993081) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6602064) q[1];
sx q[1];
rz(-1.8567137) q[1];
sx q[1];
rz(-1.6260765) q[1];
rz(-pi) q[2];
rz(-0.015720856) q[3];
sx q[3];
rz(-2.0831046) q[3];
sx q[3];
rz(-2.6920464) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.9989495) q[2];
sx q[2];
rz(-2.8012186) q[2];
sx q[2];
rz(-1.1738698) q[2];
rz(-3.0657892) q[3];
sx q[3];
rz(-1.9971763) q[3];
sx q[3];
rz(-0.092806667) q[3];
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
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3409815) q[0];
sx q[0];
rz(-1.0656463) q[0];
sx q[0];
rz(3.0766292) q[0];
rz(-2.5669572) q[1];
sx q[1];
rz(-0.42962933) q[1];
sx q[1];
rz(1.2423135) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.500538) q[0];
sx q[0];
rz(-0.28028742) q[0];
sx q[0];
rz(2.6602402) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3556446) q[2];
sx q[2];
rz(-0.64385507) q[2];
sx q[2];
rz(-1.7807963) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.2989267) q[1];
sx q[1];
rz(-1.0416404) q[1];
sx q[1];
rz(-1.3859205) q[1];
rz(-pi) q[2];
rz(-0.55450704) q[3];
sx q[3];
rz(-0.79075659) q[3];
sx q[3];
rz(0.33695541) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.3339281) q[2];
sx q[2];
rz(-2.0662722) q[2];
sx q[2];
rz(0.58369613) q[2];
rz(-2.5675473) q[3];
sx q[3];
rz(-2.0160926) q[3];
sx q[3];
rz(0.13124245) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42049256) q[0];
sx q[0];
rz(-2.2476966) q[0];
sx q[0];
rz(-0.72845355) q[0];
rz(-1.4942253) q[1];
sx q[1];
rz(-2.7431226) q[1];
sx q[1];
rz(1.0167936) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25123888) q[0];
sx q[0];
rz(-1.4892502) q[0];
sx q[0];
rz(-1.606705) q[0];
x q[1];
rz(1.0810034) q[2];
sx q[2];
rz(-1.0989597) q[2];
sx q[2];
rz(-2.8420198) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.2336725) q[1];
sx q[1];
rz(-0.83041149) q[1];
sx q[1];
rz(2.9552712) q[1];
x q[2];
rz(-3.0924762) q[3];
sx q[3];
rz(-1.8521063) q[3];
sx q[3];
rz(1.0055055) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.3399405) q[2];
sx q[2];
rz(-1.6093971) q[2];
sx q[2];
rz(2.0920848) q[2];
rz(0.63878757) q[3];
sx q[3];
rz(-2.5103266) q[3];
sx q[3];
rz(1.1857741) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3291572) q[0];
sx q[0];
rz(-1.8742467) q[0];
sx q[0];
rz(-1.4720434) q[0];
rz(2.4064348) q[1];
sx q[1];
rz(-0.77886326) q[1];
sx q[1];
rz(0.24681117) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.033337489) q[0];
sx q[0];
rz(-0.70972432) q[0];
sx q[0];
rz(1.8002585) q[0];
x q[1];
rz(2.7505789) q[2];
sx q[2];
rz(-2.6303929) q[2];
sx q[2];
rz(2.9875987) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.81891638) q[1];
sx q[1];
rz(-0.37441844) q[1];
sx q[1];
rz(0.14426343) q[1];
rz(-2.187192) q[3];
sx q[3];
rz(-1.8064926) q[3];
sx q[3];
rz(0.97782048) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.6639158) q[2];
sx q[2];
rz(-1.1185948) q[2];
sx q[2];
rz(1.6003312) q[2];
rz(0.70704308) q[3];
sx q[3];
rz(-1.0975081) q[3];
sx q[3];
rz(2.2533806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8108869) q[0];
sx q[0];
rz(-2.4208477) q[0];
sx q[0];
rz(1.8141618) q[0];
rz(1.5785626) q[1];
sx q[1];
rz(-0.47416082) q[1];
sx q[1];
rz(-0.24838233) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43863338) q[0];
sx q[0];
rz(-0.54134936) q[0];
sx q[0];
rz(0.066141733) q[0];
x q[1];
rz(2.9551198) q[2];
sx q[2];
rz(-1.7261793) q[2];
sx q[2];
rz(-1.1764256) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.76368139) q[1];
sx q[1];
rz(-1.3994201) q[1];
sx q[1];
rz(2.548449) q[1];
rz(1.8099269) q[3];
sx q[3];
rz(-2.2056747) q[3];
sx q[3];
rz(2.2374416) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.7053232) q[2];
sx q[2];
rz(-0.98781172) q[2];
sx q[2];
rz(2.3441337) q[2];
rz(2.752839) q[3];
sx q[3];
rz(-0.60384408) q[3];
sx q[3];
rz(-2.6388772) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0080863) q[0];
sx q[0];
rz(-3.0651423) q[0];
sx q[0];
rz(-1.7957934) q[0];
rz(1.0812409) q[1];
sx q[1];
rz(-1.2370279) q[1];
sx q[1];
rz(0.1246917) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.161675) q[0];
sx q[0];
rz(-1.7470164) q[0];
sx q[0];
rz(-1.4833223) q[0];
x q[1];
rz(0.6090392) q[2];
sx q[2];
rz(-1.8134724) q[2];
sx q[2];
rz(0.86953029) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.4253937) q[1];
sx q[1];
rz(-2.4024706) q[1];
sx q[1];
rz(-0.083934099) q[1];
rz(-pi) q[2];
rz(0.74495875) q[3];
sx q[3];
rz(-0.31674851) q[3];
sx q[3];
rz(-1.6572286) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.6283915) q[2];
sx q[2];
rz(-1.9548364) q[2];
sx q[2];
rz(-2.52264) q[2];
rz(1.0533054) q[3];
sx q[3];
rz(-0.17799938) q[3];
sx q[3];
rz(-0.9296023) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58182794) q[0];
sx q[0];
rz(-1.7511837) q[0];
sx q[0];
rz(-2.0986309) q[0];
rz(2.6783121) q[1];
sx q[1];
rz(-1.1136585) q[1];
sx q[1];
rz(1.0707062) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8541504) q[0];
sx q[0];
rz(-1.5702015) q[0];
sx q[0];
rz(-0.48405148) q[0];
rz(-2.7942065) q[2];
sx q[2];
rz(-1.4457448) q[2];
sx q[2];
rz(1.320968) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.067498265) q[1];
sx q[1];
rz(-2.3405582) q[1];
sx q[1];
rz(-1.4317516) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5911063) q[3];
sx q[3];
rz(-1.9029402) q[3];
sx q[3];
rz(-2.9796245) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.7730007) q[2];
sx q[2];
rz(-0.7395145) q[2];
sx q[2];
rz(-2.8179742) q[2];
rz(-0.98179022) q[3];
sx q[3];
rz(-0.86172813) q[3];
sx q[3];
rz(-2.0543082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6535646) q[0];
sx q[0];
rz(-1.7249148) q[0];
sx q[0];
rz(1.9352242) q[0];
rz(-1.9288829) q[1];
sx q[1];
rz(-0.85314631) q[1];
sx q[1];
rz(0.94747296) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0757785) q[0];
sx q[0];
rz(-1.4697207) q[0];
sx q[0];
rz(-1.0321192) q[0];
rz(-pi) q[1];
rz(2.410789) q[2];
sx q[2];
rz(-0.23554221) q[2];
sx q[2];
rz(-1.3256324) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.84711134) q[1];
sx q[1];
rz(-2.1335019) q[1];
sx q[1];
rz(-0.28114762) q[1];
rz(-pi) q[2];
x q[2];
rz(0.34668563) q[3];
sx q[3];
rz(-2.3035435) q[3];
sx q[3];
rz(-2.8447062) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.56132135) q[2];
sx q[2];
rz(-1.7611971) q[2];
sx q[2];
rz(-2.6718111) q[2];
rz(-1.8404768) q[3];
sx q[3];
rz(-1.4353292) q[3];
sx q[3];
rz(0.27967134) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25205055) q[0];
sx q[0];
rz(-0.38953504) q[0];
sx q[0];
rz(1.3289733) q[0];
rz(-0.7912311) q[1];
sx q[1];
rz(-0.33114094) q[1];
sx q[1];
rz(0.20283094) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3664353) q[0];
sx q[0];
rz(-0.046910722) q[0];
sx q[0];
rz(2.5527918) q[0];
rz(-pi) q[1];
rz(2.9183396) q[2];
sx q[2];
rz(-1.3666144) q[2];
sx q[2];
rz(0.60165652) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.087346615) q[1];
sx q[1];
rz(-0.11212238) q[1];
sx q[1];
rz(2.0263158) q[1];
rz(-2.8392302) q[3];
sx q[3];
rz(-2.6609169) q[3];
sx q[3];
rz(2.9722948) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.41708502) q[2];
sx q[2];
rz(-2.8736726) q[2];
sx q[2];
rz(2.1450796) q[2];
rz(-0.35342446) q[3];
sx q[3];
rz(-2.3963908) q[3];
sx q[3];
rz(-0.80741185) q[3];
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
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.080169454) q[0];
sx q[0];
rz(-0.81289476) q[0];
sx q[0];
rz(-0.18173519) q[0];
rz(-3.0985447) q[1];
sx q[1];
rz(-0.64518607) q[1];
sx q[1];
rz(-0.28082401) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8107287) q[0];
sx q[0];
rz(-0.64635902) q[0];
sx q[0];
rz(-0.75673639) q[0];
rz(2.4091987) q[2];
sx q[2];
rz(-2.8576982) q[2];
sx q[2];
rz(0.64901272) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.22051375) q[1];
sx q[1];
rz(-1.5346569) q[1];
sx q[1];
rz(1.5131348) q[1];
rz(-pi) q[2];
x q[2];
rz(0.19631581) q[3];
sx q[3];
rz(-1.2442949) q[3];
sx q[3];
rz(-1.6895837) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.8250371) q[2];
sx q[2];
rz(-1.8871769) q[2];
sx q[2];
rz(-2.5184856) q[2];
rz(-1.0021707) q[3];
sx q[3];
rz(-1.3391756) q[3];
sx q[3];
rz(0.56308693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4939209) q[0];
sx q[0];
rz(-1.5734084) q[0];
sx q[0];
rz(-1.5403803) q[0];
rz(-2.2676246) q[1];
sx q[1];
rz(-1.0653492) q[1];
sx q[1];
rz(-3.031562) q[1];
rz(-0.95332425) q[2];
sx q[2];
rz(-0.8669903) q[2];
sx q[2];
rz(0.16624761) q[2];
rz(2.7846746) q[3];
sx q[3];
rz(-1.1834984) q[3];
sx q[3];
rz(-0.055565861) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
