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
rz(-0.7402339) q[0];
sx q[0];
rz(-1.482168) q[0];
sx q[0];
rz(2.8066714) q[0];
rz(-2.6236293) q[1];
sx q[1];
rz(-2.1393175) q[1];
sx q[1];
rz(-0.60751539) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7142732) q[0];
sx q[0];
rz(-1.8022984) q[0];
sx q[0];
rz(1.0697068) q[0];
x q[1];
rz(3.0246441) q[2];
sx q[2];
rz(-1.0643355) q[2];
sx q[2];
rz(-3.1053271) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.40558896) q[1];
sx q[1];
rz(-1.0936671) q[1];
sx q[1];
rz(1.3911584) q[1];
x q[2];
rz(-0.3729975) q[3];
sx q[3];
rz(-1.1524876) q[3];
sx q[3];
rz(-1.4137063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.5241549) q[2];
sx q[2];
rz(-1.9257156) q[2];
sx q[2];
rz(-2.4784135) q[2];
rz(3.0607306) q[3];
sx q[3];
rz(-2.9359449) q[3];
sx q[3];
rz(-1.1801571) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2422159) q[0];
sx q[0];
rz(-1.3726534) q[0];
sx q[0];
rz(-2.2826165) q[0];
rz(-1.8513177) q[1];
sx q[1];
rz(-1.6764418) q[1];
sx q[1];
rz(1.4069517) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0075571) q[0];
sx q[0];
rz(-1.3430183) q[0];
sx q[0];
rz(-2.0640949) q[0];
rz(-pi) q[1];
x q[1];
rz(1.67703) q[2];
sx q[2];
rz(-2.8076545) q[2];
sx q[2];
rz(-0.33423697) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.3000112) q[1];
sx q[1];
rz(-1.661013) q[1];
sx q[1];
rz(2.0848227) q[1];
x q[2];
rz(-1.9858667) q[3];
sx q[3];
rz(-1.0336116) q[3];
sx q[3];
rz(1.9946757) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0591639) q[2];
sx q[2];
rz(-2.6291206) q[2];
sx q[2];
rz(1.814369) q[2];
rz(2.0969157) q[3];
sx q[3];
rz(-0.71377126) q[3];
sx q[3];
rz(2.7248342) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5663719) q[0];
sx q[0];
rz(-2.2307668) q[0];
sx q[0];
rz(2.7161993) q[0];
rz(1.7644024) q[1];
sx q[1];
rz(-1.6433989) q[1];
sx q[1];
rz(1.4345217) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4664513) q[0];
sx q[0];
rz(-1.688862) q[0];
sx q[0];
rz(2.4982014) q[0];
x q[1];
rz(-0.70830958) q[2];
sx q[2];
rz(-1.5454486) q[2];
sx q[2];
rz(-0.011653221) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.7075478) q[1];
sx q[1];
rz(-0.92869379) q[1];
sx q[1];
rz(-2.0944067) q[1];
x q[2];
rz(-1.1528963) q[3];
sx q[3];
rz(-2.4800081) q[3];
sx q[3];
rz(-2.083287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.7003358) q[2];
sx q[2];
rz(-1.8852899) q[2];
sx q[2];
rz(-1.2909935) q[2];
rz(1.357632) q[3];
sx q[3];
rz(-1.5057526) q[3];
sx q[3];
rz(-1.4972081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62035471) q[0];
sx q[0];
rz(-2.1339895) q[0];
sx q[0];
rz(0.77019101) q[0];
rz(-2.1417446) q[1];
sx q[1];
rz(-0.60037535) q[1];
sx q[1];
rz(1.8005449) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5504549) q[0];
sx q[0];
rz(-1.3584175) q[0];
sx q[0];
rz(-0.55340931) q[0];
rz(-pi) q[1];
rz(-2.6821939) q[2];
sx q[2];
rz(-2.2041577) q[2];
sx q[2];
rz(0.49672302) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.8902258) q[1];
sx q[1];
rz(-0.40491762) q[1];
sx q[1];
rz(0.7271073) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.099319066) q[3];
sx q[3];
rz(-1.2867974) q[3];
sx q[3];
rz(-1.8530396) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.3406713) q[2];
sx q[2];
rz(-1.069671) q[2];
sx q[2];
rz(-2.3731903) q[2];
rz(3.0411804) q[3];
sx q[3];
rz(-1.5407591) q[3];
sx q[3];
rz(1.8628619) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.95555821) q[0];
sx q[0];
rz(-0.30534196) q[0];
sx q[0];
rz(-2.9146063) q[0];
rz(1.7598033) q[1];
sx q[1];
rz(-0.58265668) q[1];
sx q[1];
rz(-1.3963799) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5433301) q[0];
sx q[0];
rz(-2.6285183) q[0];
sx q[0];
rz(2.3725879) q[0];
rz(-0.41047217) q[2];
sx q[2];
rz(-1.7232401) q[2];
sx q[2];
rz(-1.6688011) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.3815617) q[1];
sx q[1];
rz(-2.2027822) q[1];
sx q[1];
rz(1.4820335) q[1];
x q[2];
rz(-1.5377858) q[3];
sx q[3];
rz(-0.74652687) q[3];
sx q[3];
rz(-1.7611662) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.22623006) q[2];
sx q[2];
rz(-1.9910944) q[2];
sx q[2];
rz(-0.076233141) q[2];
rz(1.3876312) q[3];
sx q[3];
rz(-0.97532719) q[3];
sx q[3];
rz(1.7414198) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6607894) q[0];
sx q[0];
rz(-1.5810409) q[0];
sx q[0];
rz(-1.8060818) q[0];
rz(-0.90323365) q[1];
sx q[1];
rz(-1.8025554) q[1];
sx q[1];
rz(-0.18641557) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6064998) q[0];
sx q[0];
rz(-0.71487521) q[0];
sx q[0];
rz(0.8755279) q[0];
rz(0.57045649) q[2];
sx q[2];
rz(-1.4028616) q[2];
sx q[2];
rz(2.8569375) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.0707782) q[1];
sx q[1];
rz(-0.94453393) q[1];
sx q[1];
rz(0.29919821) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3535054) q[3];
sx q[3];
rz(-0.65207802) q[3];
sx q[3];
rz(1.8502082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.67941252) q[2];
sx q[2];
rz(-1.1932411) q[2];
sx q[2];
rz(1.3235486) q[2];
rz(1.1421674) q[3];
sx q[3];
rz(-0.78502941) q[3];
sx q[3];
rz(-0.89046684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(2.6776176) q[0];
sx q[0];
rz(-0.1828201) q[0];
sx q[0];
rz(0.63419813) q[0];
rz(1.0448666) q[1];
sx q[1];
rz(-2.2506782) q[1];
sx q[1];
rz(2.1814836) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7564108) q[0];
sx q[0];
rz(-2.7419081) q[0];
sx q[0];
rz(-1.6413641) q[0];
rz(-pi) q[1];
rz(3.1158218) q[2];
sx q[2];
rz(-2.8404245) q[2];
sx q[2];
rz(-2.0370551) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.83852488) q[1];
sx q[1];
rz(-2.6160598) q[1];
sx q[1];
rz(1.5938894) q[1];
x q[2];
rz(-2.3784901) q[3];
sx q[3];
rz(-2.1402485) q[3];
sx q[3];
rz(-0.72009898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1977957) q[2];
sx q[2];
rz(-2.4891977) q[2];
sx q[2];
rz(1.5459527) q[2];
rz(-1.4173077) q[3];
sx q[3];
rz(-2.3167819) q[3];
sx q[3];
rz(1.6212757) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5328131) q[0];
sx q[0];
rz(-0.19278917) q[0];
sx q[0];
rz(0.97484318) q[0];
rz(-0.10803647) q[1];
sx q[1];
rz(-1.8861176) q[1];
sx q[1];
rz(-1.9727762) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2377396) q[0];
sx q[0];
rz(-1.0971591) q[0];
sx q[0];
rz(-2.6698378) q[0];
rz(-pi) q[1];
rz(1.4270876e-05) q[2];
sx q[2];
rz(-2.37974) q[2];
sx q[2];
rz(0.95214168) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.16202422) q[1];
sx q[1];
rz(-1.8037533) q[1];
sx q[1];
rz(3.0515025) q[1];
x q[2];
rz(3.1049012) q[3];
sx q[3];
rz(-0.49202575) q[3];
sx q[3];
rz(-1.2855315) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(3.03269) q[2];
sx q[2];
rz(-2.2638075) q[2];
sx q[2];
rz(2.4857793) q[2];
rz(0.10410318) q[3];
sx q[3];
rz(-1.7963573) q[3];
sx q[3];
rz(2.9534705) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8769237) q[0];
sx q[0];
rz(-1.3035362) q[0];
sx q[0];
rz(-1.6498097) q[0];
rz(2.1022294) q[1];
sx q[1];
rz(-0.84016687) q[1];
sx q[1];
rz(-2.6731491) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.428498) q[0];
sx q[0];
rz(-1.5665496) q[0];
sx q[0];
rz(-0.18761401) q[0];
x q[1];
rz(-2.6985964) q[2];
sx q[2];
rz(-1.3123242) q[2];
sx q[2];
rz(1.4696095) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.8711618) q[1];
sx q[1];
rz(-0.79177815) q[1];
sx q[1];
rz(1.7986078) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3571068) q[3];
sx q[3];
rz(-2.753559) q[3];
sx q[3];
rz(1.7540384) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.09482065) q[2];
sx q[2];
rz(-1.5609317) q[2];
sx q[2];
rz(-2.2136733) q[2];
rz(1.8038484) q[3];
sx q[3];
rz(-2.6313621) q[3];
sx q[3];
rz(-1.6524338) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
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
rz(0.054166404) q[0];
sx q[0];
rz(-1.0324284) q[0];
sx q[0];
rz(1.077865) q[0];
rz(2.7742591) q[1];
sx q[1];
rz(-1.8108188) q[1];
sx q[1];
rz(-1.8574538) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76332247) q[0];
sx q[0];
rz(-2.5978932) q[0];
sx q[0];
rz(1.0683573) q[0];
rz(-3.025829) q[2];
sx q[2];
rz(-0.9204922) q[2];
sx q[2];
rz(2.8011326) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.2312647) q[1];
sx q[1];
rz(-2.961425) q[1];
sx q[1];
rz(2.641201) q[1];
rz(-pi) q[2];
rz(-2.4068042) q[3];
sx q[3];
rz(-1.4126724) q[3];
sx q[3];
rz(1.2722335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.0014570634) q[2];
sx q[2];
rz(-0.45150253) q[2];
sx q[2];
rz(-0.4293116) q[2];
rz(-2.9863827) q[3];
sx q[3];
rz(-0.25777543) q[3];
sx q[3];
rz(1.3252873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
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
rz(-2.8975288) q[0];
sx q[0];
rz(-1.5499935) q[0];
sx q[0];
rz(1.5503379) q[0];
rz(0.86391972) q[1];
sx q[1];
rz(-2.7680631) q[1];
sx q[1];
rz(1.6815129) q[1];
rz(2.7511394) q[2];
sx q[2];
rz(-0.57885546) q[2];
sx q[2];
rz(-0.59138966) q[2];
rz(1.1170837) q[3];
sx q[3];
rz(-2.4632005) q[3];
sx q[3];
rz(1.7076422) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
