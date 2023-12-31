OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.5326795) q[0];
sx q[0];
rz(-2.764954) q[0];
sx q[0];
rz(-0.11178804) q[0];
rz(1.6821661) q[1];
sx q[1];
rz(4.7987727) q[1];
sx q[1];
rz(6.12943) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4610447) q[0];
sx q[0];
rz(-1.8616315) q[0];
sx q[0];
rz(2.6936147) q[0];
rz(-3.0646677) q[2];
sx q[2];
rz(-1.811532) q[2];
sx q[2];
rz(-1.7488637) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.0220713) q[1];
sx q[1];
rz(-1.540456) q[1];
sx q[1];
rz(1.8896709) q[1];
x q[2];
rz(-0.93011232) q[3];
sx q[3];
rz(-2.6112587) q[3];
sx q[3];
rz(-1.9213284) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.443632) q[2];
sx q[2];
rz(-1.7093095) q[2];
sx q[2];
rz(-1.4367746) q[2];
rz(-2.4076961) q[3];
sx q[3];
rz(-1.5926444) q[3];
sx q[3];
rz(2.6255887) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2519418) q[0];
sx q[0];
rz(-1.9152859) q[0];
sx q[0];
rz(0.92457986) q[0];
rz(-0.997116) q[1];
sx q[1];
rz(-0.50874248) q[1];
sx q[1];
rz(-1.3234214) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.130587) q[0];
sx q[0];
rz(-2.4936094) q[0];
sx q[0];
rz(0.76052163) q[0];
rz(-pi) q[1];
rz(1.5263444) q[2];
sx q[2];
rz(-1.2215081) q[2];
sx q[2];
rz(-2.8291707) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.70222774) q[1];
sx q[1];
rz(-0.80184466) q[1];
sx q[1];
rz(-2.761809) q[1];
rz(0.99625846) q[3];
sx q[3];
rz(-1.901445) q[3];
sx q[3];
rz(-1.656267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.9374342) q[2];
sx q[2];
rz(-1.5896475) q[2];
sx q[2];
rz(2.3834174) q[2];
rz(-2.5126863) q[3];
sx q[3];
rz(-0.40142504) q[3];
sx q[3];
rz(1.988407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73308289) q[0];
sx q[0];
rz(-1.874431) q[0];
sx q[0];
rz(2.4531903) q[0];
rz(-0.06772659) q[1];
sx q[1];
rz(-1.7522782) q[1];
sx q[1];
rz(-0.53007954) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5085174) q[0];
sx q[0];
rz(-1.8607664) q[0];
sx q[0];
rz(0.11147186) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7128503) q[2];
sx q[2];
rz(-1.0849761) q[2];
sx q[2];
rz(-0.57810099) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-3.0878151) q[1];
sx q[1];
rz(-2.3936845) q[1];
sx q[1];
rz(1.0653711) q[1];
rz(1.9534555) q[3];
sx q[3];
rz(-2.8842162) q[3];
sx q[3];
rz(0.60636683) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.80785859) q[2];
sx q[2];
rz(-0.01161751) q[2];
sx q[2];
rz(2.2401436) q[2];
rz(2.3060913) q[3];
sx q[3];
rz(-1.52799) q[3];
sx q[3];
rz(-1.8301331) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19105844) q[0];
sx q[0];
rz(-1.5427417) q[0];
sx q[0];
rz(2.3572671) q[0];
rz(-0.061231881) q[1];
sx q[1];
rz(-0.71413723) q[1];
sx q[1];
rz(3.004946) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.559583) q[0];
sx q[0];
rz(-1.0916545) q[0];
sx q[0];
rz(2.9942306) q[0];
rz(-pi) q[1];
rz(2.6776671) q[2];
sx q[2];
rz(-1.8564965) q[2];
sx q[2];
rz(-2.6505016) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.0285443) q[1];
sx q[1];
rz(-1.5946769) q[1];
sx q[1];
rz(-1.6391812) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.09282077) q[3];
sx q[3];
rz(-0.99675677) q[3];
sx q[3];
rz(2.1294347) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.0121997) q[2];
sx q[2];
rz(-0.94038525) q[2];
sx q[2];
rz(-0.56048918) q[2];
rz(-0.012332049) q[3];
sx q[3];
rz(-0.90356946) q[3];
sx q[3];
rz(-2.0509317) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43301582) q[0];
sx q[0];
rz(-2.5362159) q[0];
sx q[0];
rz(2.3204455) q[0];
rz(0.87617809) q[1];
sx q[1];
rz(-0.89996243) q[1];
sx q[1];
rz(1.7339773) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1836023) q[0];
sx q[0];
rz(-1.2215658) q[0];
sx q[0];
rz(-2.7244085) q[0];
rz(-0.3048004) q[2];
sx q[2];
rz(-1.0101724) q[2];
sx q[2];
rz(1.0501109) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.7387313) q[1];
sx q[1];
rz(-2.1364473) q[1];
sx q[1];
rz(-2.9823142) q[1];
x q[2];
rz(2.9191454) q[3];
sx q[3];
rz(-1.9223833) q[3];
sx q[3];
rz(0.44140154) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.6901107) q[2];
sx q[2];
rz(-1.2146981) q[2];
sx q[2];
rz(-0.042479854) q[2];
rz(-2.5111607) q[3];
sx q[3];
rz(-0.63215956) q[3];
sx q[3];
rz(2.650034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38934389) q[0];
sx q[0];
rz(-1.2279953) q[0];
sx q[0];
rz(-2.65843) q[0];
rz(-2.0893611) q[1];
sx q[1];
rz(-1.9960884) q[1];
sx q[1];
rz(2.5767456) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42445499) q[0];
sx q[0];
rz(-1.0466482) q[0];
sx q[0];
rz(0.15630396) q[0];
rz(3.0658423) q[2];
sx q[2];
rz(-1.9905914) q[2];
sx q[2];
rz(-2.9043353) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.7397592) q[1];
sx q[1];
rz(-1.1439011) q[1];
sx q[1];
rz(3.0763807) q[1];
rz(-pi) q[2];
rz(-1.7259898) q[3];
sx q[3];
rz(-1.3421913) q[3];
sx q[3];
rz(-2.5582563) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.5715282) q[2];
sx q[2];
rz(-1.0674942) q[2];
sx q[2];
rz(-1.8035536) q[2];
rz(-1.3048874) q[3];
sx q[3];
rz(-2.0740502) q[3];
sx q[3];
rz(0.6750955) q[3];
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
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17334443) q[0];
sx q[0];
rz(-1.4302379) q[0];
sx q[0];
rz(0.54779732) q[0];
rz(-2.3563747) q[1];
sx q[1];
rz(-1.3350057) q[1];
sx q[1];
rz(-0.26842591) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3184506) q[0];
sx q[0];
rz(-2.8822691) q[0];
sx q[0];
rz(0.32487049) q[0];
rz(-pi) q[1];
rz(2.354291) q[2];
sx q[2];
rz(-2.1848218) q[2];
sx q[2];
rz(1.3348483) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.56663471) q[1];
sx q[1];
rz(-2.1556427) q[1];
sx q[1];
rz(-3.023874) q[1];
rz(-2.1106342) q[3];
sx q[3];
rz(-1.3098048) q[3];
sx q[3];
rz(0.48418448) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.61775529) q[2];
sx q[2];
rz(-2.338151) q[2];
sx q[2];
rz(0.8141554) q[2];
rz(-0.37627775) q[3];
sx q[3];
rz(-1.9777931) q[3];
sx q[3];
rz(-3.0686839) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
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
rz(2.5557264) q[0];
sx q[0];
rz(-1.4050452) q[0];
sx q[0];
rz(2.1222173) q[0];
rz(-0.85340071) q[1];
sx q[1];
rz(-1.1420206) q[1];
sx q[1];
rz(0.44874915) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1837346) q[0];
sx q[0];
rz(-2.0439889) q[0];
sx q[0];
rz(0.26259043) q[0];
rz(-pi) q[1];
rz(2.8837187) q[2];
sx q[2];
rz(-1.2784625) q[2];
sx q[2];
rz(-0.1567947) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.110072) q[1];
sx q[1];
rz(-1.3168465) q[1];
sx q[1];
rz(-0.54540789) q[1];
x q[2];
rz(2.2711146) q[3];
sx q[3];
rz(-0.92910367) q[3];
sx q[3];
rz(0.93922797) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.86924187) q[2];
sx q[2];
rz(-1.3808455) q[2];
sx q[2];
rz(-0.31420079) q[2];
rz(-0.82434404) q[3];
sx q[3];
rz(-2.6894675) q[3];
sx q[3];
rz(2.3468988) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0712414) q[0];
sx q[0];
rz(-3.0817139) q[0];
sx q[0];
rz(-1.8810133) q[0];
rz(0.64385995) q[1];
sx q[1];
rz(-1.9088129) q[1];
sx q[1];
rz(3.1226645) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3138872) q[0];
sx q[0];
rz(-0.2642309) q[0];
sx q[0];
rz(-3.0731191) q[0];
rz(-pi) q[1];
rz(0.55982121) q[2];
sx q[2];
rz(-2.9580742) q[2];
sx q[2];
rz(0.57932094) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.1259809) q[1];
sx q[1];
rz(-0.19128448) q[1];
sx q[1];
rz(-2.8730238) q[1];
x q[2];
rz(0.67846672) q[3];
sx q[3];
rz(-2.7760091) q[3];
sx q[3];
rz(1.6775223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.54946047) q[2];
sx q[2];
rz(-2.7533054) q[2];
sx q[2];
rz(0.67031676) q[2];
rz(-0.51236764) q[3];
sx q[3];
rz(-1.7497601) q[3];
sx q[3];
rz(1.7211154) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55384127) q[0];
sx q[0];
rz(-1.8974263) q[0];
sx q[0];
rz(-0.49945369) q[0];
rz(-1.5669426) q[1];
sx q[1];
rz(-0.27856871) q[1];
sx q[1];
rz(2.0589028) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64718819) q[0];
sx q[0];
rz(-1.6196961) q[0];
sx q[0];
rz(0.58736579) q[0];
rz(2.322299) q[2];
sx q[2];
rz(-1.594992) q[2];
sx q[2];
rz(2.1404612) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.6607099) q[1];
sx q[1];
rz(-2.1143267) q[1];
sx q[1];
rz(-1.5488312) q[1];
rz(-pi) q[2];
rz(-0.38090221) q[3];
sx q[3];
rz(-1.2657832) q[3];
sx q[3];
rz(-2.0088793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.3988951) q[2];
sx q[2];
rz(-1.164914) q[2];
sx q[2];
rz(2.6297074) q[2];
rz(-2.7486457) q[3];
sx q[3];
rz(-1.7397375) q[3];
sx q[3];
rz(-1.0569364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.186541) q[0];
sx q[0];
rz(-1.3165836) q[0];
sx q[0];
rz(-2.5008428) q[0];
rz(-0.74116771) q[1];
sx q[1];
rz(-2.3186431) q[1];
sx q[1];
rz(2.9021312) q[1];
rz(-0.95460931) q[2];
sx q[2];
rz(-1.90345) q[2];
sx q[2];
rz(-0.30751139) q[2];
rz(1.3988597) q[3];
sx q[3];
rz(-0.83172432) q[3];
sx q[3];
rz(-1.5749501) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
