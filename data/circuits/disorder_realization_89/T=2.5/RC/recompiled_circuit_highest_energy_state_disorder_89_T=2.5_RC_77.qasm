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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3008904) q[0];
sx q[0];
rz(-1.7127174) q[0];
sx q[0];
rz(3.0034756) q[0];
rz(-pi) q[1];
rz(1.5643372) q[2];
sx q[2];
rz(-1.1273555) q[2];
sx q[2];
rz(2.9270767) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.2121316) q[1];
sx q[1];
rz(-0.60570215) q[1];
sx q[1];
rz(-0.49537658) q[1];
rz(1.3985083) q[3];
sx q[3];
rz(-0.43753657) q[3];
sx q[3];
rz(-3.1222536) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.67285377) q[2];
sx q[2];
rz(-1.3287975) q[2];
sx q[2];
rz(1.9358181) q[2];
rz(2.2852211) q[3];
sx q[3];
rz(-2.2090293) q[3];
sx q[3];
rz(2.2642853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67742753) q[0];
sx q[0];
rz(-0.39746118) q[0];
sx q[0];
rz(-2.933266) q[0];
rz(-1.9147929) q[1];
sx q[1];
rz(-2.0715163) q[1];
sx q[1];
rz(-0.56455451) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78075942) q[0];
sx q[0];
rz(-1.5281023) q[0];
sx q[0];
rz(-1.4111931) q[0];
rz(-pi) q[1];
rz(2.0656086) q[2];
sx q[2];
rz(-2.2391265) q[2];
sx q[2];
rz(0.47023103) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.7617088) q[1];
sx q[1];
rz(-1.392388) q[1];
sx q[1];
rz(-2.4316923) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6030231) q[3];
sx q[3];
rz(-2.6015266) q[3];
sx q[3];
rz(-2.9750254) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.63689256) q[2];
sx q[2];
rz(-1.2496639) q[2];
sx q[2];
rz(-1.6553817) q[2];
rz(-0.49556035) q[3];
sx q[3];
rz(-1.7335745) q[3];
sx q[3];
rz(2.1347617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(0.14690742) q[0];
sx q[0];
rz(-2.4201396) q[0];
sx q[0];
rz(-1.7578693) q[0];
rz(-1.9231298) q[1];
sx q[1];
rz(-2.0530901) q[1];
sx q[1];
rz(0.4247492) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.527194) q[0];
sx q[0];
rz(-1.5081811) q[0];
sx q[0];
rz(1.6024496) q[0];
rz(1.3176383) q[2];
sx q[2];
rz(-1.8033529) q[2];
sx q[2];
rz(-0.81040224) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.9806165) q[1];
sx q[1];
rz(-2.7572639) q[1];
sx q[1];
rz(0.20341441) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5307755) q[3];
sx q[3];
rz(-0.95453018) q[3];
sx q[3];
rz(1.5933526) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.72033) q[2];
sx q[2];
rz(-2.0739372) q[2];
sx q[2];
rz(3.1171411) q[2];
rz(-1.9931741) q[3];
sx q[3];
rz(-2.0423856) q[3];
sx q[3];
rz(1.0816157) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9352683) q[0];
sx q[0];
rz(-0.22980389) q[0];
sx q[0];
rz(2.9277053) q[0];
rz(0.81854406) q[1];
sx q[1];
rz(-1.7844424) q[1];
sx q[1];
rz(0.68665409) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.548107) q[0];
sx q[0];
rz(-1.0229551) q[0];
sx q[0];
rz(2.0440897) q[0];
rz(-pi) q[1];
rz(-2.1193723) q[2];
sx q[2];
rz(-1.3201664) q[2];
sx q[2];
rz(-0.20285367) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.80685341) q[1];
sx q[1];
rz(-0.76419965) q[1];
sx q[1];
rz(-0.13607009) q[1];
x q[2];
rz(-2.2591001) q[3];
sx q[3];
rz(-0.63418856) q[3];
sx q[3];
rz(2.2435482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.116918) q[2];
sx q[2];
rz(-1.7322098) q[2];
sx q[2];
rz(-2.4521258) q[2];
rz(1.3245448) q[3];
sx q[3];
rz(-1.2169415) q[3];
sx q[3];
rz(2.569258) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6784742) q[0];
sx q[0];
rz(-0.093570396) q[0];
sx q[0];
rz(2.1043188) q[0];
rz(-2.4689238) q[1];
sx q[1];
rz(-1.172352) q[1];
sx q[1];
rz(2.345828) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0507495) q[0];
sx q[0];
rz(-1.2552869) q[0];
sx q[0];
rz(2.1054224) q[0];
rz(-pi) q[1];
rz(-0.99822171) q[2];
sx q[2];
rz(-1.8279982) q[2];
sx q[2];
rz(2.5679979) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.9059674) q[1];
sx q[1];
rz(-0.86938953) q[1];
sx q[1];
rz(2.7011306) q[1];
rz(0.66821258) q[3];
sx q[3];
rz(-1.3283806) q[3];
sx q[3];
rz(-0.53430623) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.1611288) q[2];
sx q[2];
rz(-1.0995862) q[2];
sx q[2];
rz(-2.3058092) q[2];
rz(-0.25660822) q[3];
sx q[3];
rz(-2.0308688) q[3];
sx q[3];
rz(-0.49032828) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68818727) q[0];
sx q[0];
rz(-1.9904933) q[0];
sx q[0];
rz(1.9292462) q[0];
rz(-1.0351828) q[1];
sx q[1];
rz(-1.6500708) q[1];
sx q[1];
rz(1.0119247) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8261004) q[0];
sx q[0];
rz(-1.6796711) q[0];
sx q[0];
rz(-1.7427765) q[0];
rz(-pi) q[1];
rz(1.6056608) q[2];
sx q[2];
rz(-2.8375585) q[2];
sx q[2];
rz(-2.3140098) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.1862667) q[1];
sx q[1];
rz(-1.1623993) q[1];
sx q[1];
rz(-0.46640654) q[1];
rz(-pi) q[2];
rz(2.5402734) q[3];
sx q[3];
rz(-2.6868002) q[3];
sx q[3];
rz(-3.0704344) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6637491) q[2];
sx q[2];
rz(-0.19890824) q[2];
sx q[2];
rz(-0.94856962) q[2];
rz(-2.6712724) q[3];
sx q[3];
rz(-2.0430653) q[3];
sx q[3];
rz(-1.2547913) q[3];
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
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54742852) q[0];
sx q[0];
rz(-2.1305888) q[0];
sx q[0];
rz(-1.8773361) q[0];
rz(2.3612379) q[1];
sx q[1];
rz(-2.5006313) q[1];
sx q[1];
rz(-0.19187127) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.980968) q[0];
sx q[0];
rz(-1.1066529) q[0];
sx q[0];
rz(1.1726517) q[0];
rz(-pi) q[1];
rz(-1.1926417) q[2];
sx q[2];
rz(-1.9072235) q[2];
sx q[2];
rz(1.9371197) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.33436066) q[1];
sx q[1];
rz(-1.8340115) q[1];
sx q[1];
rz(-0.46137793) q[1];
rz(-pi) q[2];
rz(2.7924813) q[3];
sx q[3];
rz(-0.99945949) q[3];
sx q[3];
rz(-0.1196564) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.4912305) q[2];
sx q[2];
rz(-1.0543084) q[2];
sx q[2];
rz(-0.80797705) q[2];
rz(-0.84336495) q[3];
sx q[3];
rz(-1.9873514) q[3];
sx q[3];
rz(1.556506) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(1.4845487) q[0];
sx q[0];
rz(-3.1175685) q[0];
sx q[0];
rz(-2.3727544) q[0];
rz(-0.27278849) q[1];
sx q[1];
rz(-1.5130679) q[1];
sx q[1];
rz(1.7441033) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5834806) q[0];
sx q[0];
rz(-1.2654103) q[0];
sx q[0];
rz(0.38174133) q[0];
rz(-pi) q[1];
rz(-2.0968123) q[2];
sx q[2];
rz(-2.5616841) q[2];
sx q[2];
rz(-1.2666262) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.8500415) q[1];
sx q[1];
rz(-1.631389) q[1];
sx q[1];
rz(-1.687595) q[1];
rz(0.52941771) q[3];
sx q[3];
rz(-2.4000958) q[3];
sx q[3];
rz(-1.3138258) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.2046795) q[2];
sx q[2];
rz(-2.0759373) q[2];
sx q[2];
rz(1.5020471) q[2];
rz(-1.3192734) q[3];
sx q[3];
rz(-0.84984142) q[3];
sx q[3];
rz(-1.5892861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69916344) q[0];
sx q[0];
rz(-0.84446877) q[0];
sx q[0];
rz(2.1635639) q[0];
rz(2.6507822) q[1];
sx q[1];
rz(-1.699479) q[1];
sx q[1];
rz(-1.4960272) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7154244) q[0];
sx q[0];
rz(-1.4472571) q[0];
sx q[0];
rz(0.23793313) q[0];
x q[1];
rz(-2.7649058) q[2];
sx q[2];
rz(-0.44038195) q[2];
sx q[2];
rz(2.1747766) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.5514976) q[1];
sx q[1];
rz(-1.2480772) q[1];
sx q[1];
rz(1.3151874) q[1];
x q[2];
rz(0.53197143) q[3];
sx q[3];
rz(-1.5051748) q[3];
sx q[3];
rz(-2.402879) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.28828037) q[2];
sx q[2];
rz(-2.2744982) q[2];
sx q[2];
rz(-1.8193998) q[2];
rz(-0.36618048) q[3];
sx q[3];
rz(-1.7145232) q[3];
sx q[3];
rz(2.0100994) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2243097) q[0];
sx q[0];
rz(-1.6095105) q[0];
sx q[0];
rz(3.068058) q[0];
rz(-1.3167943) q[1];
sx q[1];
rz(-1.2578745) q[1];
sx q[1];
rz(-1.8427461) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4816581) q[0];
sx q[0];
rz(-1.3941688) q[0];
sx q[0];
rz(-0.28089471) q[0];
rz(-pi) q[1];
rz(-1.3234947) q[2];
sx q[2];
rz(-0.50116936) q[2];
sx q[2];
rz(0.71401087) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9292752) q[1];
sx q[1];
rz(-1.3354509) q[1];
sx q[1];
rz(-1.1849684) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4019626) q[3];
sx q[3];
rz(-1.6508023) q[3];
sx q[3];
rz(0.39749872) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.1765882) q[2];
sx q[2];
rz(-2.1829288) q[2];
sx q[2];
rz(0.57903543) q[2];
rz(-1.6252919) q[3];
sx q[3];
rz(-2.5987891) q[3];
sx q[3];
rz(-1.496284) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2753006) q[0];
sx q[0];
rz(-1.0466812) q[0];
sx q[0];
rz(-0.064591797) q[0];
rz(0.91118377) q[1];
sx q[1];
rz(-1.8235527) q[1];
sx q[1];
rz(1.8240737) q[1];
rz(0.0077132465) q[2];
sx q[2];
rz(-2.7174502) q[2];
sx q[2];
rz(2.9677718) q[2];
rz(0.28476247) q[3];
sx q[3];
rz(-0.97281572) q[3];
sx q[3];
rz(0.90481467) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
