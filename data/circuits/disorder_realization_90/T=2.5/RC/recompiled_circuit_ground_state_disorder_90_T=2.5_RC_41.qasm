OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.5034135) q[0];
sx q[0];
rz(-1.2794275) q[0];
sx q[0];
rz(-2.3655565) q[0];
rz(-2.4322721) q[1];
sx q[1];
rz(-1.5172989) q[1];
sx q[1];
rz(-2.5425743) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.020928362) q[0];
sx q[0];
rz(-1.0987704) q[0];
sx q[0];
rz(-0.13735227) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8593596) q[2];
sx q[2];
rz(-2.2405409) q[2];
sx q[2];
rz(1.9053659) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.4645965) q[1];
sx q[1];
rz(-1.076816) q[1];
sx q[1];
rz(-1.9736273) q[1];
rz(-0.78035155) q[3];
sx q[3];
rz(-1.8858741) q[3];
sx q[3];
rz(1.9596069) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.1385931) q[2];
sx q[2];
rz(-2.004576) q[2];
sx q[2];
rz(-2.9453759) q[2];
rz(-1.044322) q[3];
sx q[3];
rz(-2.315867) q[3];
sx q[3];
rz(2.830982) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0394548) q[0];
sx q[0];
rz(-2.4882443) q[0];
sx q[0];
rz(2.2838604) q[0];
rz(-0.54025447) q[1];
sx q[1];
rz(-2.0876355) q[1];
sx q[1];
rz(2.3589755) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58843553) q[0];
sx q[0];
rz(-0.62214506) q[0];
sx q[0];
rz(1.4216656) q[0];
rz(-2.906053) q[2];
sx q[2];
rz(-0.96220926) q[2];
sx q[2];
rz(-0.71266251) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.5021584) q[1];
sx q[1];
rz(-2.5089988) q[1];
sx q[1];
rz(2.8117287) q[1];
x q[2];
rz(-0.55763839) q[3];
sx q[3];
rz(-2.1774315) q[3];
sx q[3];
rz(-2.5733657) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.1648078) q[2];
sx q[2];
rz(-0.79443496) q[2];
sx q[2];
rz(0.59107333) q[2];
rz(1.0278541) q[3];
sx q[3];
rz(-0.63575345) q[3];
sx q[3];
rz(-0.13636057) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22981055) q[0];
sx q[0];
rz(-0.52017838) q[0];
sx q[0];
rz(-1.0562563) q[0];
rz(-2.1504869) q[1];
sx q[1];
rz(-1.7678363) q[1];
sx q[1];
rz(0.80449218) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1209542) q[0];
sx q[0];
rz(-2.197663) q[0];
sx q[0];
rz(1.2406209) q[0];
rz(-pi) q[1];
rz(-0.94203888) q[2];
sx q[2];
rz(-1.4938746) q[2];
sx q[2];
rz(-0.49401894) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.3884405) q[1];
sx q[1];
rz(-0.33666753) q[1];
sx q[1];
rz(-2.3364288) q[1];
x q[2];
rz(0.64739703) q[3];
sx q[3];
rz(-1.7248099) q[3];
sx q[3];
rz(1.2274418) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.450401) q[2];
sx q[2];
rz(-2.2914026) q[2];
sx q[2];
rz(-2.5737393) q[2];
rz(1.9495226) q[3];
sx q[3];
rz(-2.6046533) q[3];
sx q[3];
rz(-0.15538628) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6464226) q[0];
sx q[0];
rz(-1.7429054) q[0];
sx q[0];
rz(-1.1871185) q[0];
rz(-2.8861956) q[1];
sx q[1];
rz(-1.4030158) q[1];
sx q[1];
rz(2.752221) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0856649) q[0];
sx q[0];
rz(-1.4367668) q[0];
sx q[0];
rz(-2.9197951) q[0];
x q[1];
rz(1.6045447) q[2];
sx q[2];
rz(-0.72740388) q[2];
sx q[2];
rz(2.6836065) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.8592718) q[1];
sx q[1];
rz(-0.92073694) q[1];
sx q[1];
rz(-0.23843022) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7674753) q[3];
sx q[3];
rz(-1.4405319) q[3];
sx q[3];
rz(-1.2127339) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.4046459) q[2];
sx q[2];
rz(-0.32361042) q[2];
sx q[2];
rz(1.3673937) q[2];
rz(0.5395475) q[3];
sx q[3];
rz(-1.5893785) q[3];
sx q[3];
rz(1.3246983) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0747727) q[0];
sx q[0];
rz(-0.19817752) q[0];
sx q[0];
rz(-0.5603801) q[0];
rz(-2.9226774) q[1];
sx q[1];
rz(-1.0812662) q[1];
sx q[1];
rz(0.17094368) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6534916) q[0];
sx q[0];
rz(-0.91918901) q[0];
sx q[0];
rz(2.1733858) q[0];
rz(-1.7124699) q[2];
sx q[2];
rz(-1.3686841) q[2];
sx q[2];
rz(2.0703482) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.39300181) q[1];
sx q[1];
rz(-0.95470631) q[1];
sx q[1];
rz(2.2115117) q[1];
rz(-1.4724949) q[3];
sx q[3];
rz(-0.62255961) q[3];
sx q[3];
rz(-0.05427256) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.35974744) q[2];
sx q[2];
rz(-2.013194) q[2];
sx q[2];
rz(-2.183059) q[2];
rz(-0.16252276) q[3];
sx q[3];
rz(-2.2992117) q[3];
sx q[3];
rz(-1.5490279) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.573134) q[0];
sx q[0];
rz(-3.077226) q[0];
sx q[0];
rz(-2.3934613) q[0];
rz(2.023078) q[1];
sx q[1];
rz(-2.7225814) q[1];
sx q[1];
rz(0.31707877) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7359809) q[0];
sx q[0];
rz(-1.6389585) q[0];
sx q[0];
rz(0.24644417) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.19642475) q[2];
sx q[2];
rz(-1.6187917) q[2];
sx q[2];
rz(-0.53507198) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.78817716) q[1];
sx q[1];
rz(-1.4995575) q[1];
sx q[1];
rz(-2.337238) q[1];
rz(-3.1012332) q[3];
sx q[3];
rz(-0.77733126) q[3];
sx q[3];
rz(-3.0632927) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.3743484) q[2];
sx q[2];
rz(-2.2129462) q[2];
sx q[2];
rz(1.7283776) q[2];
rz(-3.0610541) q[3];
sx q[3];
rz(-1.7929411) q[3];
sx q[3];
rz(-0.22182375) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(0.11248511) q[0];
sx q[0];
rz(-1.8928098) q[0];
sx q[0];
rz(-1.2741733) q[0];
rz(1.796465) q[1];
sx q[1];
rz(-2.3990264) q[1];
sx q[1];
rz(1.8003731) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8250834) q[0];
sx q[0];
rz(-2.1232455) q[0];
sx q[0];
rz(0.98487206) q[0];
x q[1];
rz(-2.6244782) q[2];
sx q[2];
rz(-0.54143751) q[2];
sx q[2];
rz(-1.1351897) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.53645027) q[1];
sx q[1];
rz(-1.4172232) q[1];
sx q[1];
rz(-1.3773514) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2209211) q[3];
sx q[3];
rz(-1.6643644) q[3];
sx q[3];
rz(2.2232995) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.207927) q[2];
sx q[2];
rz(-1.1873446) q[2];
sx q[2];
rz(-2.6742317) q[2];
rz(0.76006877) q[3];
sx q[3];
rz(-1.4773388) q[3];
sx q[3];
rz(1.2490341) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6811328) q[0];
sx q[0];
rz(-2.12119) q[0];
sx q[0];
rz(-1.4469752) q[0];
rz(0.32577062) q[1];
sx q[1];
rz(-0.21369801) q[1];
sx q[1];
rz(-1.6206585) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5012799) q[0];
sx q[0];
rz(-2.5569943) q[0];
sx q[0];
rz(2.6452097) q[0];
rz(1.5704186) q[2];
sx q[2];
rz(-1.6624358) q[2];
sx q[2];
rz(0.82110345) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.28315) q[1];
sx q[1];
rz(-1.9236739) q[1];
sx q[1];
rz(-0.18650413) q[1];
rz(-1.8160024) q[3];
sx q[3];
rz(-0.71906656) q[3];
sx q[3];
rz(1.2326425) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.1763566) q[2];
sx q[2];
rz(-0.74702817) q[2];
sx q[2];
rz(-2.526324) q[2];
rz(-1.8687013) q[3];
sx q[3];
rz(-2.8764909) q[3];
sx q[3];
rz(0.42241514) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1212921) q[0];
sx q[0];
rz(-1.8676119) q[0];
sx q[0];
rz(-0.85025382) q[0];
rz(-1.1212768) q[1];
sx q[1];
rz(-1.3669776) q[1];
sx q[1];
rz(0.77176315) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92895011) q[0];
sx q[0];
rz(-1.54689) q[0];
sx q[0];
rz(2.2646409) q[0];
rz(-pi) q[1];
rz(2.9608126) q[2];
sx q[2];
rz(-1.5018936) q[2];
sx q[2];
rz(3.1071752) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.96442879) q[1];
sx q[1];
rz(-1.5156974) q[1];
sx q[1];
rz(-1.3269618) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8949299) q[3];
sx q[3];
rz(-1.3132902) q[3];
sx q[3];
rz(2.3070564) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-3.0261592) q[2];
sx q[2];
rz(-1.5449646) q[2];
sx q[2];
rz(2.3010632) q[2];
rz(-3.0606411) q[3];
sx q[3];
rz(-2.5637124) q[3];
sx q[3];
rz(0.24278434) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42846546) q[0];
sx q[0];
rz(-2.9436538) q[0];
sx q[0];
rz(-1.9848829) q[0];
rz(-2.1139862) q[1];
sx q[1];
rz(-1.7637858) q[1];
sx q[1];
rz(0.81046945) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5549703) q[0];
sx q[0];
rz(-1.7867309) q[0];
sx q[0];
rz(2.3761889) q[0];
x q[1];
rz(1.2357622) q[2];
sx q[2];
rz(-2.0910237) q[2];
sx q[2];
rz(-0.65641415) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.1406882) q[1];
sx q[1];
rz(-1.8546805) q[1];
sx q[1];
rz(1.8769916) q[1];
rz(2.9083071) q[3];
sx q[3];
rz(-0.57396997) q[3];
sx q[3];
rz(-3.1149394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.79601866) q[2];
sx q[2];
rz(-1.2556262) q[2];
sx q[2];
rz(-1.6176809) q[2];
rz(1.1003305) q[3];
sx q[3];
rz(-1.9610145) q[3];
sx q[3];
rz(-2.9617917) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0236459) q[0];
sx q[0];
rz(-1.4143586) q[0];
sx q[0];
rz(2.216862) q[0];
rz(-0.042451518) q[1];
sx q[1];
rz(-2.0274542) q[1];
sx q[1];
rz(-1.7995119) q[1];
rz(1.7003851) q[2];
sx q[2];
rz(-2.7322506) q[2];
sx q[2];
rz(2.6553497) q[2];
rz(3.068559) q[3];
sx q[3];
rz(-0.45847736) q[3];
sx q[3];
rz(1.7401742) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
