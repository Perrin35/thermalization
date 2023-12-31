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
rz(3.5182313) q[0];
sx q[0];
rz(9.3129899) q[0];
rz(1.6821661) q[1];
sx q[1];
rz(4.7987727) q[1];
sx q[1];
rz(6.12943) q[1];
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
rz(-0.15969294) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.598688) q[1];
sx q[1];
rz(-0.3202657) q[1];
sx q[1];
rz(1.6673052) q[1];
x q[2];
rz(0.93011232) q[3];
sx q[3];
rz(-2.6112587) q[3];
sx q[3];
rz(1.9213284) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.6979606) q[2];
sx q[2];
rz(-1.7093095) q[2];
sx q[2];
rz(-1.704818) q[2];
rz(2.4076961) q[3];
sx q[3];
rz(-1.5926444) q[3];
sx q[3];
rz(-2.6255887) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2519418) q[0];
sx q[0];
rz(-1.2263068) q[0];
sx q[0];
rz(-0.92457986) q[0];
rz(-0.997116) q[1];
sx q[1];
rz(-2.6328502) q[1];
sx q[1];
rz(1.3234214) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2307777) q[0];
sx q[0];
rz(-1.1417023) q[0];
sx q[0];
rz(-0.50165117) q[0];
rz(1.6152482) q[2];
sx q[2];
rz(-1.9200846) q[2];
sx q[2];
rz(-2.8291707) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.70222774) q[1];
sx q[1];
rz(-0.80184466) q[1];
sx q[1];
rz(2.761809) q[1];
x q[2];
rz(-0.38815659) q[3];
sx q[3];
rz(-1.0309439) q[3];
sx q[3];
rz(2.848958) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.9374342) q[2];
sx q[2];
rz(-1.5519451) q[2];
sx q[2];
rz(2.3834174) q[2];
rz(-2.5126863) q[3];
sx q[3];
rz(-2.7401676) q[3];
sx q[3];
rz(1.1531856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.73308289) q[0];
sx q[0];
rz(-1.2671616) q[0];
sx q[0];
rz(-0.68840233) q[0];
rz(-3.0738661) q[1];
sx q[1];
rz(-1.3893145) q[1];
sx q[1];
rz(-0.53007954) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90572689) q[0];
sx q[0];
rz(-1.6775963) q[0];
sx q[0];
rz(1.279116) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0968139) q[2];
sx q[2];
rz(-1.194343) q[2];
sx q[2];
rz(1.9386171) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.5487919) q[1];
sx q[1];
rz(-0.93344102) q[1];
sx q[1];
rz(-2.7194276) q[1];
rz(-0.097966627) q[3];
sx q[3];
rz(-1.8091822) q[3];
sx q[3];
rz(1.0007678) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.3337341) q[2];
sx q[2];
rz(-3.1299751) q[2];
sx q[2];
rz(-0.90144908) q[2];
rz(2.3060913) q[3];
sx q[3];
rz(-1.52799) q[3];
sx q[3];
rz(-1.8301331) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19105844) q[0];
sx q[0];
rz(-1.598851) q[0];
sx q[0];
rz(0.78432551) q[0];
rz(-3.0803608) q[1];
sx q[1];
rz(-2.4274554) q[1];
sx q[1];
rz(3.004946) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58200968) q[0];
sx q[0];
rz(-2.0499381) q[0];
sx q[0];
rz(2.9942306) q[0];
rz(-2.5606974) q[2];
sx q[2];
rz(-0.53933203) q[2];
sx q[2];
rz(2.5748594) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.0285443) q[1];
sx q[1];
rz(-1.5469157) q[1];
sx q[1];
rz(1.6391812) q[1];
x q[2];
rz(1.7131545) q[3];
sx q[3];
rz(-2.5609303) q[3];
sx q[3];
rz(-2.2992087) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.1293929) q[2];
sx q[2];
rz(-0.94038525) q[2];
sx q[2];
rz(-2.5811035) q[2];
rz(-3.1292606) q[3];
sx q[3];
rz(-2.2380232) q[3];
sx q[3];
rz(-2.0509317) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(2.7085768) q[0];
sx q[0];
rz(-2.5362159) q[0];
sx q[0];
rz(-0.82114712) q[0];
rz(2.2654146) q[1];
sx q[1];
rz(-0.89996243) q[1];
sx q[1];
rz(1.4076153) q[1];
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
rz(2.1528835) q[2];
sx q[2];
rz(-1.313813) q[2];
sx q[2];
rz(-0.35494057) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.40286139) q[1];
sx q[1];
rz(-1.0051454) q[1];
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
rz(0.63043198) q[3];
sx q[3];
rz(-2.5094331) q[3];
sx q[3];
rz(0.49155864) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38934389) q[0];
sx q[0];
rz(-1.9135973) q[0];
sx q[0];
rz(-2.65843) q[0];
rz(-1.0522316) q[1];
sx q[1];
rz(-1.9960884) q[1];
sx q[1];
rz(0.56484708) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7171377) q[0];
sx q[0];
rz(-2.0949445) q[0];
sx q[0];
rz(0.15630396) q[0];
x q[1];
rz(-1.4028366) q[2];
sx q[2];
rz(-0.42617455) q[2];
sx q[2];
rz(3.0884398) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.7397592) q[1];
sx q[1];
rz(-1.1439011) q[1];
sx q[1];
rz(-3.0763807) q[1];
rz(-pi) q[2];
rz(-1.4156028) q[3];
sx q[3];
rz(-1.7994013) q[3];
sx q[3];
rz(-2.5582563) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.57006449) q[2];
sx q[2];
rz(-1.0674942) q[2];
sx q[2];
rz(-1.338039) q[2];
rz(-1.3048874) q[3];
sx q[3];
rz(-1.0675425) q[3];
sx q[3];
rz(2.4664972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9682482) q[0];
sx q[0];
rz(-1.7113547) q[0];
sx q[0];
rz(-2.5937953) q[0];
rz(-2.3563747) q[1];
sx q[1];
rz(-1.806587) q[1];
sx q[1];
rz(-2.8731667) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3184506) q[0];
sx q[0];
rz(-0.25932352) q[0];
sx q[0];
rz(0.32487049) q[0];
rz(-0.7873017) q[2];
sx q[2];
rz(-0.95677081) q[2];
sx q[2];
rz(1.8067443) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.35560289) q[1];
sx q[1];
rz(-0.59521788) q[1];
sx q[1];
rz(1.7463513) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0309585) q[3];
sx q[3];
rz(-1.8317878) q[3];
sx q[3];
rz(2.6574082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.5238374) q[2];
sx q[2];
rz(-2.338151) q[2];
sx q[2];
rz(-0.8141554) q[2];
rz(0.37627775) q[3];
sx q[3];
rz(-1.9777931) q[3];
sx q[3];
rz(-0.072908727) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58586621) q[0];
sx q[0];
rz(-1.7365475) q[0];
sx q[0];
rz(2.1222173) q[0];
rz(2.2881919) q[1];
sx q[1];
rz(-1.9995721) q[1];
sx q[1];
rz(2.6928435) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42496029) q[0];
sx q[0];
rz(-2.6052931) q[0];
sx q[0];
rz(2.0400356) q[0];
rz(-pi) q[1];
rz(-0.25787392) q[2];
sx q[2];
rz(-1.2784625) q[2];
sx q[2];
rz(2.984798) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.9948506) q[1];
sx q[1];
rz(-0.59616201) q[1];
sx q[1];
rz(-2.6776828) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3677424) q[3];
sx q[3];
rz(-1.0283096) q[3];
sx q[3];
rz(2.0427996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.2723508) q[2];
sx q[2];
rz(-1.7607471) q[2];
sx q[2];
rz(0.31420079) q[2];
rz(2.3172486) q[3];
sx q[3];
rz(-2.6894675) q[3];
sx q[3];
rz(2.3468988) q[3];
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
x q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.070351275) q[0];
sx q[0];
rz(-3.0817139) q[0];
sx q[0];
rz(-1.2605793) q[0];
rz(-0.64385995) q[1];
sx q[1];
rz(-1.2327797) q[1];
sx q[1];
rz(-0.018928122) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3138872) q[0];
sx q[0];
rz(-2.8773617) q[0];
sx q[0];
rz(0.068473579) q[0];
rz(-pi) q[1];
rz(1.47255) q[2];
sx q[2];
rz(-1.7260523) q[2];
sx q[2];
rz(-1.9948024) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.2889274) q[1];
sx q[1];
rz(-1.3864494) q[1];
sx q[1];
rz(-1.6221371) q[1];
rz(-pi) q[2];
x q[2];
rz(0.28963611) q[3];
sx q[3];
rz(-1.3445065) q[3];
sx q[3];
rz(-0.7520523) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.5921322) q[2];
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
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55384127) q[0];
sx q[0];
rz(-1.8974263) q[0];
sx q[0];
rz(-2.642139) q[0];
rz(-1.5669426) q[1];
sx q[1];
rz(-0.27856871) q[1];
sx q[1];
rz(-1.0826899) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4944045) q[0];
sx q[0];
rz(-1.6196961) q[0];
sx q[0];
rz(-0.58736579) q[0];
rz(-pi) q[1];
rz(0.033109025) q[2];
sx q[2];
rz(-0.81956714) q[2];
sx q[2];
rz(-2.5493252) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.078552695) q[1];
sx q[1];
rz(-1.5895956) q[1];
sx q[1];
rz(-0.54363721) q[1];
rz(-pi) q[2];
rz(-2.7606904) q[3];
sx q[3];
rz(-1.2657832) q[3];
sx q[3];
rz(-1.1327133) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.7426976) q[2];
sx q[2];
rz(-1.9766786) q[2];
sx q[2];
rz(2.6297074) q[2];
rz(2.7486457) q[3];
sx q[3];
rz(-1.7397375) q[3];
sx q[3];
rz(-2.0846562) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95505161) q[0];
sx q[0];
rz(-1.3165836) q[0];
sx q[0];
rz(-2.5008428) q[0];
rz(-2.4004249) q[1];
sx q[1];
rz(-0.82294958) q[1];
sx q[1];
rz(-0.23946147) q[1];
rz(-2.1869833) q[2];
sx q[2];
rz(-1.2381427) q[2];
sx q[2];
rz(2.8340813) q[2];
rz(-2.9560271) q[3];
sx q[3];
rz(-0.75510988) q[3];
sx q[3];
rz(1.818944) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
