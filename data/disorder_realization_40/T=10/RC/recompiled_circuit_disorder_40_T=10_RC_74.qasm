OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.6089132) q[0];
sx q[0];
rz(-0.37663868) q[0];
sx q[0];
rz(-3.0298046) q[0];
rz(-1.4594266) q[1];
sx q[1];
rz(-1.6571801) q[1];
sx q[1];
rz(-2.9878374) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.789334) q[0];
sx q[0];
rz(-2.6129122) q[0];
sx q[0];
rz(2.5369011) q[0];
rz(-pi) q[1];
rz(-1.8122187) q[2];
sx q[2];
rz(-1.6454988) q[2];
sx q[2];
rz(2.9818997) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.54290463) q[1];
sx q[1];
rz(-2.821327) q[1];
sx q[1];
rz(-1.4742875) q[1];
x q[2];
rz(2.2114803) q[3];
sx q[3];
rz(-2.6112587) q[3];
sx q[3];
rz(-1.9213284) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.6979606) q[2];
sx q[2];
rz(-1.4322832) q[2];
sx q[2];
rz(1.4367746) q[2];
rz(0.73389655) q[3];
sx q[3];
rz(-1.5489483) q[3];
sx q[3];
rz(-2.6255887) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.88965082) q[0];
sx q[0];
rz(-1.9152859) q[0];
sx q[0];
rz(-0.92457986) q[0];
rz(-2.1444767) q[1];
sx q[1];
rz(-2.6328502) q[1];
sx q[1];
rz(-1.3234214) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.01100563) q[0];
sx q[0];
rz(-0.6479833) q[0];
sx q[0];
rz(0.76052163) q[0];
rz(-2.7919865) q[2];
sx q[2];
rz(-1.5290302) q[2];
sx q[2];
rz(-1.2431527) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.18113187) q[1];
sx q[1];
rz(-0.84003969) q[1];
sx q[1];
rz(1.204927) q[1];
rz(2.7534361) q[3];
sx q[3];
rz(-2.1106488) q[3];
sx q[3];
rz(-2.848958) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.20415846) q[2];
sx q[2];
rz(-1.5519451) q[2];
sx q[2];
rz(2.3834174) q[2];
rz(-0.6289064) q[3];
sx q[3];
rz(-0.40142504) q[3];
sx q[3];
rz(-1.988407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73308289) q[0];
sx q[0];
rz(-1.2671616) q[0];
sx q[0];
rz(0.68840233) q[0];
rz(3.0738661) q[1];
sx q[1];
rz(-1.7522782) q[1];
sx q[1];
rz(2.6115131) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.006223) q[0];
sx q[0];
rz(-2.831499) q[0];
sx q[0];
rz(-1.2139411) q[0];
x q[1];
rz(2.0968139) q[2];
sx q[2];
rz(-1.9472497) q[2];
sx q[2];
rz(1.9386171) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.9024132) q[1];
sx q[1];
rz(-1.2352408) q[1];
sx q[1];
rz(-2.2526342) q[1];
rz(-pi) q[2];
rz(-1.1881371) q[3];
sx q[3];
rz(-0.2573765) q[3];
sx q[3];
rz(-0.60636683) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.80785859) q[2];
sx q[2];
rz(-3.1299751) q[2];
sx q[2];
rz(0.90144908) q[2];
rz(-2.3060913) q[3];
sx q[3];
rz(-1.52799) q[3];
sx q[3];
rz(1.8301331) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9505342) q[0];
sx q[0];
rz(-1.5427417) q[0];
sx q[0];
rz(-2.3572671) q[0];
rz(0.061231881) q[1];
sx q[1];
rz(-0.71413723) q[1];
sx q[1];
rz(-3.004946) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89350677) q[0];
sx q[0];
rz(-2.6419905) q[0];
sx q[0];
rz(1.2953555) q[0];
rz(1.253445) q[2];
sx q[2];
rz(-1.1270521) q[2];
sx q[2];
rz(-1.2197942) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.0285443) q[1];
sx q[1];
rz(-1.5946769) q[1];
sx q[1];
rz(-1.5024115) q[1];
x q[2];
rz(1.7131545) q[3];
sx q[3];
rz(-2.5609303) q[3];
sx q[3];
rz(-2.2992087) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.0121997) q[2];
sx q[2];
rz(-2.2012074) q[2];
sx q[2];
rz(-0.56048918) q[2];
rz(0.012332049) q[3];
sx q[3];
rz(-0.90356946) q[3];
sx q[3];
rz(-1.0906609) q[3];
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
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7085768) q[0];
sx q[0];
rz(-2.5362159) q[0];
sx q[0];
rz(-0.82114712) q[0];
rz(-0.87617809) q[1];
sx q[1];
rz(-0.89996243) q[1];
sx q[1];
rz(-1.7339773) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1836023) q[0];
sx q[0];
rz(-1.2215658) q[0];
sx q[0];
rz(2.7244085) q[0];
rz(-2.0166964) q[2];
sx q[2];
rz(-0.6302399) q[2];
sx q[2];
rz(-1.5843887) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.059543) q[1];
sx q[1];
rz(-1.4364916) q[1];
sx q[1];
rz(2.1422269) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9191454) q[3];
sx q[3];
rz(-1.2192093) q[3];
sx q[3];
rz(2.7001911) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.451482) q[2];
sx q[2];
rz(-1.2146981) q[2];
sx q[2];
rz(-0.042479854) q[2];
rz(2.5111607) q[3];
sx q[3];
rz(-2.5094331) q[3];
sx q[3];
rz(2.650034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7522488) q[0];
sx q[0];
rz(-1.2279953) q[0];
sx q[0];
rz(2.65843) q[0];
rz(-2.0893611) q[1];
sx q[1];
rz(-1.1455043) q[1];
sx q[1];
rz(0.56484708) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42445499) q[0];
sx q[0];
rz(-2.0949445) q[0];
sx q[0];
rz(-0.15630396) q[0];
rz(1.1499314) q[2];
sx q[2];
rz(-1.5016342) q[2];
sx q[2];
rz(1.7771306) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.5833252) q[1];
sx q[1];
rz(-0.43154432) q[1];
sx q[1];
rz(1.4285018) q[1];
rz(-pi) q[2];
rz(-1.7259898) q[3];
sx q[3];
rz(-1.3421913) q[3];
sx q[3];
rz(-2.5582563) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.5715282) q[2];
sx q[2];
rz(-2.0740985) q[2];
sx q[2];
rz(1.8035536) q[2];
rz(1.8367052) q[3];
sx q[3];
rz(-2.0740502) q[3];
sx q[3];
rz(-2.4664972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9682482) q[0];
sx q[0];
rz(-1.7113547) q[0];
sx q[0];
rz(-0.54779732) q[0];
rz(-2.3563747) q[1];
sx q[1];
rz(-1.3350057) q[1];
sx q[1];
rz(2.8731667) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9831532) q[0];
sx q[0];
rz(-1.8162677) q[0];
sx q[0];
rz(-1.4863187) q[0];
rz(0.78598778) q[2];
sx q[2];
rz(-2.18835) q[2];
sx q[2];
rz(0.76030375) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.5749579) q[1];
sx q[1];
rz(-2.1556427) q[1];
sx q[1];
rz(3.023874) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0915756) q[3];
sx q[3];
rz(-2.5476533) q[3];
sx q[3];
rz(-1.493243) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5238374) q[2];
sx q[2];
rz(-0.80344168) q[2];
sx q[2];
rz(-0.8141554) q[2];
rz(-0.37627775) q[3];
sx q[3];
rz(-1.1637996) q[3];
sx q[3];
rz(3.0686839) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58586621) q[0];
sx q[0];
rz(-1.7365475) q[0];
sx q[0];
rz(1.0193753) q[0];
rz(0.85340071) q[1];
sx q[1];
rz(-1.1420206) q[1];
sx q[1];
rz(-0.44874915) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9578581) q[0];
sx q[0];
rz(-2.0439889) q[0];
sx q[0];
rz(0.26259043) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2690527) q[2];
sx q[2];
rz(-1.3240959) q[2];
sx q[2];
rz(-1.6517284) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.1467421) q[1];
sx q[1];
rz(-0.59616201) q[1];
sx q[1];
rz(-0.46390987) q[1];
rz(2.3677424) q[3];
sx q[3];
rz(-2.1132831) q[3];
sx q[3];
rz(-2.0427996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.86924187) q[2];
sx q[2];
rz(-1.7607471) q[2];
sx q[2];
rz(-2.8273919) q[2];
rz(-2.3172486) q[3];
sx q[3];
rz(-2.6894675) q[3];
sx q[3];
rz(-2.3468988) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0712414) q[0];
sx q[0];
rz(-0.059878778) q[0];
sx q[0];
rz(1.8810133) q[0];
rz(0.64385995) q[1];
sx q[1];
rz(-1.2327797) q[1];
sx q[1];
rz(0.018928122) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32301329) q[0];
sx q[0];
rz(-1.5529263) q[0];
sx q[0];
rz(-2.8779526) q[0];
rz(2.9855965) q[2];
sx q[2];
rz(-1.6678572) q[2];
sx q[2];
rz(-2.7023466) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.1259809) q[1];
sx q[1];
rz(-0.19128448) q[1];
sx q[1];
rz(0.26856883) q[1];
rz(-pi) q[2];
rz(-2.8519565) q[3];
sx q[3];
rz(-1.3445065) q[3];
sx q[3];
rz(2.3895404) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.54946047) q[2];
sx q[2];
rz(-0.38828725) q[2];
sx q[2];
rz(-2.4712759) q[2];
rz(0.51236764) q[3];
sx q[3];
rz(-1.3918326) q[3];
sx q[3];
rz(-1.4204773) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55384127) q[0];
sx q[0];
rz(-1.8974263) q[0];
sx q[0];
rz(-2.642139) q[0];
rz(1.5746501) q[1];
sx q[1];
rz(-2.8630239) q[1];
sx q[1];
rz(-2.0589028) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4944045) q[0];
sx q[0];
rz(-1.5218966) q[0];
sx q[0];
rz(-0.58736579) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.81929368) q[2];
sx q[2];
rz(-1.594992) q[2];
sx q[2];
rz(-1.0011315) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.078552695) q[1];
sx q[1];
rz(-1.551997) q[1];
sx q[1];
rz(-2.5979554) q[1];
rz(1.8977676) q[3];
sx q[3];
rz(-1.2083112) q[3];
sx q[3];
rz(2.8231951) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.7426976) q[2];
sx q[2];
rz(-1.164914) q[2];
sx q[2];
rz(2.6297074) q[2];
rz(-2.7486457) q[3];
sx q[3];
rz(-1.7397375) q[3];
sx q[3];
rz(2.0846562) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.186541) q[0];
sx q[0];
rz(-1.825009) q[0];
sx q[0];
rz(0.64074989) q[0];
rz(2.4004249) q[1];
sx q[1];
rz(-2.3186431) q[1];
sx q[1];
rz(2.9021312) q[1];
rz(-1.0319866) q[2];
sx q[2];
rz(-2.451755) q[2];
sx q[2];
rz(-2.3103726) q[2];
rz(-1.742733) q[3];
sx q[3];
rz(-0.83172432) q[3];
sx q[3];
rz(-1.5749501) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];