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
rz(0.97857082) q[0];
sx q[0];
rz(-0.11712722) q[0];
sx q[0];
rz(0.22476619) q[0];
rz(2.6449142) q[1];
sx q[1];
rz(4.7159046) q[1];
sx q[1];
rz(11.386303) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4311572) q[0];
sx q[0];
rz(-1.4340766) q[0];
sx q[0];
rz(-1.7140634) q[0];
rz(-pi) q[1];
rz(2.6981437) q[2];
sx q[2];
rz(-1.5649619) q[2];
sx q[2];
rz(1.782541) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.22332033) q[1];
sx q[1];
rz(-1.8448571) q[1];
sx q[1];
rz(2.5943701) q[1];
rz(-pi) q[2];
rz(-1.7430844) q[3];
sx q[3];
rz(-2.7040561) q[3];
sx q[3];
rz(-0.019339081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.4687389) q[2];
sx q[2];
rz(-1.8127952) q[2];
sx q[2];
rz(1.2057745) q[2];
rz(-2.2852211) q[3];
sx q[3];
rz(-2.2090293) q[3];
sx q[3];
rz(-2.2642853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(2.4641651) q[0];
sx q[0];
rz(-0.39746118) q[0];
sx q[0];
rz(0.20832668) q[0];
rz(1.9147929) q[1];
sx q[1];
rz(-2.0715163) q[1];
sx q[1];
rz(0.56455451) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3608332) q[0];
sx q[0];
rz(-1.6134904) q[0];
sx q[0];
rz(-1.4111931) q[0];
rz(-pi) q[1];
rz(0.73123572) q[2];
sx q[2];
rz(-1.1888767) q[2];
sx q[2];
rz(-0.77788355) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.98726749) q[1];
sx q[1];
rz(-2.4134153) q[1];
sx q[1];
rz(-2.8716692) q[1];
x q[2];
rz(-1.0309593) q[3];
sx q[3];
rz(-1.5542277) q[3];
sx q[3];
rz(1.709721) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.63689256) q[2];
sx q[2];
rz(-1.8919287) q[2];
sx q[2];
rz(1.6553817) q[2];
rz(2.6460323) q[3];
sx q[3];
rz(-1.4080181) q[3];
sx q[3];
rz(-2.1347617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
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
rz(0.14690742) q[0];
sx q[0];
rz(-2.4201396) q[0];
sx q[0];
rz(1.3837234) q[0];
rz(1.9231298) q[1];
sx q[1];
rz(-2.0530901) q[1];
sx q[1];
rz(-0.4247492) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61439862) q[0];
sx q[0];
rz(-1.6334115) q[0];
sx q[0];
rz(-1.539143) q[0];
rz(-pi) q[1];
rz(0.23992691) q[2];
sx q[2];
rz(-1.3245965) q[2];
sx q[2];
rz(2.4407516) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.9207117) q[1];
sx q[1];
rz(-1.6466116) q[1];
sx q[1];
rz(-2.7644509) q[1];
rz(-pi) q[2];
rz(-2.5249486) q[3];
sx q[3];
rz(-1.5381406) q[3];
sx q[3];
rz(0.045696229) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.72033) q[2];
sx q[2];
rz(-1.0676554) q[2];
sx q[2];
rz(0.024451582) q[2];
rz(1.1484185) q[3];
sx q[3];
rz(-1.0992071) q[3];
sx q[3];
rz(-1.0816157) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20632437) q[0];
sx q[0];
rz(-0.22980389) q[0];
sx q[0];
rz(0.21388737) q[0];
rz(0.81854406) q[1];
sx q[1];
rz(-1.7844424) q[1];
sx q[1];
rz(-2.4549386) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71664602) q[0];
sx q[0];
rz(-1.970463) q[0];
sx q[0];
rz(2.5406688) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8501007) q[2];
sx q[2];
rz(-1.0412058) q[2];
sx q[2];
rz(1.6232217) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.5221646) q[1];
sx q[1];
rz(-2.3261737) q[1];
sx q[1];
rz(-1.700089) q[1];
rz(-pi) q[2];
rz(2.0874168) q[3];
sx q[3];
rz(-1.1849019) q[3];
sx q[3];
rz(-1.2579045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.116918) q[2];
sx q[2];
rz(-1.4093829) q[2];
sx q[2];
rz(0.68946687) q[2];
rz(1.8170478) q[3];
sx q[3];
rz(-1.2169415) q[3];
sx q[3];
rz(0.57233468) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4631184) q[0];
sx q[0];
rz(-3.0480223) q[0];
sx q[0];
rz(-1.0372739) q[0];
rz(2.4689238) q[1];
sx q[1];
rz(-1.172352) q[1];
sx q[1];
rz(0.79576463) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0507495) q[0];
sx q[0];
rz(-1.8863057) q[0];
sx q[0];
rz(2.1054224) q[0];
rz(-pi) q[1];
rz(0.3032844) q[2];
sx q[2];
rz(-1.0192843) q[2];
sx q[2];
rz(-1.9818652) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.2356253) q[1];
sx q[1];
rz(-0.86938953) q[1];
sx q[1];
rz(-2.7011306) q[1];
x q[2];
rz(0.37973399) q[3];
sx q[3];
rz(-0.70443166) q[3];
sx q[3];
rz(-1.80988) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.98046389) q[2];
sx q[2];
rz(-2.0420065) q[2];
sx q[2];
rz(-2.3058092) q[2];
rz(-0.25660822) q[3];
sx q[3];
rz(-2.0308688) q[3];
sx q[3];
rz(2.6512644) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68818727) q[0];
sx q[0];
rz(-1.1510993) q[0];
sx q[0];
rz(-1.9292462) q[0];
rz(-1.0351828) q[1];
sx q[1];
rz(-1.4915219) q[1];
sx q[1];
rz(-1.0119247) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8261004) q[0];
sx q[0];
rz(-1.6796711) q[0];
sx q[0];
rz(-1.7427765) q[0];
rz(-pi) q[1];
x q[1];
rz(3.1306562) q[2];
sx q[2];
rz(-1.8746398) q[2];
sx q[2];
rz(-2.2774709) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.95532596) q[1];
sx q[1];
rz(-1.9791934) q[1];
sx q[1];
rz(0.46640654) q[1];
rz(-pi) q[2];
rz(2.5402734) q[3];
sx q[3];
rz(-2.6868002) q[3];
sx q[3];
rz(-3.0704344) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.6637491) q[2];
sx q[2];
rz(-2.9426844) q[2];
sx q[2];
rz(-2.193023) q[2];
rz(-2.6712724) q[3];
sx q[3];
rz(-1.0985273) q[3];
sx q[3];
rz(-1.8868014) q[3];
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
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5941641) q[0];
sx q[0];
rz(-2.1305888) q[0];
sx q[0];
rz(-1.2642566) q[0];
rz(-0.7803548) q[1];
sx q[1];
rz(-2.5006313) q[1];
sx q[1];
rz(-0.19187127) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3656552) q[0];
sx q[0];
rz(-1.2167551) q[0];
sx q[0];
rz(2.6440622) q[0];
x q[1];
rz(1.1926417) q[2];
sx q[2];
rz(-1.9072235) q[2];
sx q[2];
rz(1.2044729) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.75429186) q[1];
sx q[1];
rz(-2.6151492) q[1];
sx q[1];
rz(2.597288) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0818008) q[3];
sx q[3];
rz(-2.4823175) q[3];
sx q[3];
rz(2.6687572) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.65036217) q[2];
sx q[2];
rz(-2.0872842) q[2];
sx q[2];
rz(-0.80797705) q[2];
rz(-2.2982277) q[3];
sx q[3];
rz(-1.1542412) q[3];
sx q[3];
rz(-1.5850867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4845487) q[0];
sx q[0];
rz(-0.024024155) q[0];
sx q[0];
rz(0.76883823) q[0];
rz(-2.8688042) q[1];
sx q[1];
rz(-1.6285248) q[1];
sx q[1];
rz(1.7441033) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0087939) q[0];
sx q[0];
rz(-1.2075612) q[0];
sx q[0];
rz(1.2433267) q[0];
rz(-1.0447803) q[2];
sx q[2];
rz(-2.5616841) q[2];
sx q[2];
rz(1.2666262) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.7556909) q[1];
sx q[1];
rz(-3.0100757) q[1];
sx q[1];
rz(2.0507858) q[1];
rz(2.0040183) q[3];
sx q[3];
rz(-2.1931291) q[3];
sx q[3];
rz(-2.4985) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.93691319) q[2];
sx q[2];
rz(-2.0759373) q[2];
sx q[2];
rz(-1.6395456) q[2];
rz(1.3192734) q[3];
sx q[3];
rz(-2.2917512) q[3];
sx q[3];
rz(1.5523065) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69916344) q[0];
sx q[0];
rz(-2.2971239) q[0];
sx q[0];
rz(0.97802877) q[0];
rz(2.6507822) q[1];
sx q[1];
rz(-1.4421137) q[1];
sx q[1];
rz(-1.6455654) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.967088) q[0];
sx q[0];
rz(-1.8068815) q[0];
sx q[0];
rz(-1.4437136) q[0];
rz(2.7285887) q[2];
sx q[2];
rz(-1.7282515) q[2];
sx q[2];
rz(0.26034376) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.2435808) q[1];
sx q[1];
rz(-1.812938) q[1];
sx q[1];
rz(-0.33278521) q[1];
rz(-pi) q[2];
x q[2];
rz(0.53197143) q[3];
sx q[3];
rz(-1.6364179) q[3];
sx q[3];
rz(-0.73871368) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.8533123) q[2];
sx q[2];
rz(-2.2744982) q[2];
sx q[2];
rz(1.8193998) q[2];
rz(-2.7754122) q[3];
sx q[3];
rz(-1.4270695) q[3];
sx q[3];
rz(2.0100994) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
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
rz(-2.2243097) q[0];
sx q[0];
rz(-1.5320822) q[0];
sx q[0];
rz(-0.073534615) q[0];
rz(-1.8247983) q[1];
sx q[1];
rz(-1.8837181) q[1];
sx q[1];
rz(-1.8427461) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5054437) q[0];
sx q[0];
rz(-2.8110286) q[0];
sx q[0];
rz(0.57204582) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3234947) q[2];
sx q[2];
rz(-2.6404233) q[2];
sx q[2];
rz(0.71401087) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.9292752) q[1];
sx q[1];
rz(-1.3354509) q[1];
sx q[1];
rz(-1.9566243) q[1];
rz(2.4019626) q[3];
sx q[3];
rz(-1.4907903) q[3];
sx q[3];
rz(0.39749872) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.9650044) q[2];
sx q[2];
rz(-0.95866385) q[2];
sx q[2];
rz(-2.5625572) q[2];
rz(-1.6252919) q[3];
sx q[3];
rz(-2.5987891) q[3];
sx q[3];
rz(1.6453086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8662921) q[0];
sx q[0];
rz(-2.0949114) q[0];
sx q[0];
rz(3.0770009) q[0];
rz(0.91118377) q[1];
sx q[1];
rz(-1.8235527) q[1];
sx q[1];
rz(1.8240737) q[1];
rz(-0.0077132465) q[2];
sx q[2];
rz(-0.42414244) q[2];
sx q[2];
rz(-0.17382081) q[2];
rz(-2.1880423) q[3];
sx q[3];
rz(-1.336477) q[3];
sx q[3];
rz(-0.82930641) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
