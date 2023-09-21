OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.2620579) q[0];
sx q[0];
rz(7.6927778) q[0];
sx q[0];
rz(11.132244) q[0];
rz(0.6342451) q[1];
sx q[1];
rz(6.8847818) q[1];
sx q[1];
rz(9.8431982) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2375862) q[0];
sx q[0];
rz(-1.8338086) q[0];
sx q[0];
rz(2.0531274) q[0];
x q[1];
rz(-0.11169545) q[2];
sx q[2];
rz(-1.4354435) q[2];
sx q[2];
rz(0.40458194) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.1474255) q[1];
sx q[1];
rz(-2.6174183) q[1];
sx q[1];
rz(2.0736573) q[1];
rz(-pi) q[2];
rz(-2.9207346) q[3];
sx q[3];
rz(-1.1410032) q[3];
sx q[3];
rz(2.5151099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.3216386) q[2];
sx q[2];
rz(-1.7724089) q[2];
sx q[2];
rz(-0.83797541) q[2];
rz(2.6485802) q[3];
sx q[3];
rz(-2.8686782) q[3];
sx q[3];
rz(-3.0626007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
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
rz(0.74719602) q[0];
sx q[0];
rz(-2.4173739) q[0];
sx q[0];
rz(1.2778506) q[0];
rz(0.17678075) q[1];
sx q[1];
rz(-1.3143833) q[1];
sx q[1];
rz(-0.4321672) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5837001) q[0];
sx q[0];
rz(-1.0520792) q[0];
sx q[0];
rz(1.8899263) q[0];
rz(-1.1439267) q[2];
sx q[2];
rz(-1.9383213) q[2];
sx q[2];
rz(1.2183684) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.6501112) q[1];
sx q[1];
rz(-2.0108129) q[1];
sx q[1];
rz(3.0122709) q[1];
rz(-pi) q[2];
rz(-1.6099986) q[3];
sx q[3];
rz(-0.56534846) q[3];
sx q[3];
rz(0.70985868) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.5923578) q[2];
sx q[2];
rz(-1.2640415) q[2];
sx q[2];
rz(2.6548927) q[2];
rz(-1.7633847) q[3];
sx q[3];
rz(-1.2599726) q[3];
sx q[3];
rz(2.6087705) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.10087) q[0];
sx q[0];
rz(-0.7464872) q[0];
sx q[0];
rz(2.7242463) q[0];
rz(1.4886645) q[1];
sx q[1];
rz(-2.5960943) q[1];
sx q[1];
rz(-2.6352077) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9837512) q[0];
sx q[0];
rz(-2.7451773) q[0];
sx q[0];
rz(1.0979963) q[0];
rz(2.8767013) q[2];
sx q[2];
rz(-2.3893917) q[2];
sx q[2];
rz(1.0571935) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.80458927) q[1];
sx q[1];
rz(-2.5924006) q[1];
sx q[1];
rz(-2.5227491) q[1];
rz(1.1261602) q[3];
sx q[3];
rz(-0.4053084) q[3];
sx q[3];
rz(2.0733548) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.69904077) q[2];
sx q[2];
rz(-2.6802345) q[2];
sx q[2];
rz(-0.59147269) q[2];
rz(2.5555723) q[3];
sx q[3];
rz(-1.932671) q[3];
sx q[3];
rz(-1.4311786) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
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
rz(-1.1699003) q[0];
sx q[0];
rz(-0.62830347) q[0];
sx q[0];
rz(-2.3024094) q[0];
rz(-3.1160141) q[1];
sx q[1];
rz(-2.4459116) q[1];
sx q[1];
rz(1.5930088) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.87985932) q[0];
sx q[0];
rz(-1.121959) q[0];
sx q[0];
rz(-0.10106048) q[0];
rz(0.84393878) q[2];
sx q[2];
rz(-2.2288449) q[2];
sx q[2];
rz(2.7968) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.38494426) q[1];
sx q[1];
rz(-1.1886016) q[1];
sx q[1];
rz(-0.76809831) q[1];
rz(-pi) q[2];
rz(1.2876835) q[3];
sx q[3];
rz(-1.1554171) q[3];
sx q[3];
rz(2.0302041) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.8362391) q[2];
sx q[2];
rz(-2.2576015) q[2];
sx q[2];
rz(3.0419066) q[2];
rz(0.95885197) q[3];
sx q[3];
rz(-1.3189664) q[3];
sx q[3];
rz(1.8166186) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56617671) q[0];
sx q[0];
rz(-1.7594936) q[0];
sx q[0];
rz(-0.25594041) q[0];
rz(0.4610962) q[1];
sx q[1];
rz(-1.0436811) q[1];
sx q[1];
rz(2.3815313) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3292023) q[0];
sx q[0];
rz(-2.985552) q[0];
sx q[0];
rz(0.72326707) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0779482) q[2];
sx q[2];
rz(-1.9846989) q[2];
sx q[2];
rz(2.6457583) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.45436817) q[1];
sx q[1];
rz(-1.4441274) q[1];
sx q[1];
rz(1.441799) q[1];
rz(-pi) q[2];
x q[2];
rz(0.36565904) q[3];
sx q[3];
rz(-0.95698157) q[3];
sx q[3];
rz(2.1706276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.5710859) q[2];
sx q[2];
rz(-1.3723624) q[2];
sx q[2];
rz(2.467353) q[2];
rz(2.9267866) q[3];
sx q[3];
rz(-0.45682296) q[3];
sx q[3];
rz(-3.1242483) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6102585) q[0];
sx q[0];
rz(-1.6711618) q[0];
sx q[0];
rz(-1.9624788) q[0];
rz(-2.9367661) q[1];
sx q[1];
rz(-2.3463459) q[1];
sx q[1];
rz(2.0746322) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2444006) q[0];
sx q[0];
rz(-3.1163437) q[0];
sx q[0];
rz(0.61141725) q[0];
rz(-pi) q[1];
rz(0.42491575) q[2];
sx q[2];
rz(-1.9746466) q[2];
sx q[2];
rz(-1.0645107) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.9251717) q[1];
sx q[1];
rz(-0.74837084) q[1];
sx q[1];
rz(1.640635) q[1];
x q[2];
rz(-0.79640572) q[3];
sx q[3];
rz(-1.685433) q[3];
sx q[3];
rz(-0.20533268) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.431488) q[2];
sx q[2];
rz(-1.8523214) q[2];
sx q[2];
rz(-2.5816494) q[2];
rz(-0.7263178) q[3];
sx q[3];
rz(-2.8328219) q[3];
sx q[3];
rz(0.30549756) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85754919) q[0];
sx q[0];
rz(-2.5950268) q[0];
sx q[0];
rz(1.7204826) q[0];
rz(-2.9395318) q[1];
sx q[1];
rz(-1.4338564) q[1];
sx q[1];
rz(-0.85817671) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9265384) q[0];
sx q[0];
rz(-2.9460322) q[0];
sx q[0];
rz(-2.2901448) q[0];
rz(-pi) q[1];
x q[1];
rz(1.467534) q[2];
sx q[2];
rz(-2.729136) q[2];
sx q[2];
rz(0.39706424) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.8691751) q[1];
sx q[1];
rz(-0.032748001) q[1];
sx q[1];
rz(-2.6564024) q[1];
rz(-1.9167561) q[3];
sx q[3];
rz(-1.02709) q[3];
sx q[3];
rz(2.791415) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8537366) q[2];
sx q[2];
rz(-0.47984543) q[2];
sx q[2];
rz(1.8161592) q[2];
rz(2.251513) q[3];
sx q[3];
rz(-1.1471014) q[3];
sx q[3];
rz(2.1530698) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6778075) q[0];
sx q[0];
rz(-0.57254922) q[0];
sx q[0];
rz(2.7668787) q[0];
rz(-0.97887865) q[1];
sx q[1];
rz(-0.68190494) q[1];
sx q[1];
rz(1.3495061) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6222854) q[0];
sx q[0];
rz(-1.1847727) q[0];
sx q[0];
rz(2.4849154) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.026272341) q[2];
sx q[2];
rz(-0.71825829) q[2];
sx q[2];
rz(2.9422613) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.046557758) q[1];
sx q[1];
rz(-1.3890424) q[1];
sx q[1];
rz(-3.0912116) q[1];
rz(-pi) q[2];
rz(-1.8324864) q[3];
sx q[3];
rz(-1.040254) q[3];
sx q[3];
rz(-1.7390651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.4902041) q[2];
sx q[2];
rz(-0.75275246) q[2];
sx q[2];
rz(2.6728969) q[2];
rz(1.9474585) q[3];
sx q[3];
rz(-1.8374551) q[3];
sx q[3];
rz(1.3635925) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
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
rz(-0.63012183) q[0];
sx q[0];
rz(-2.2352495) q[0];
sx q[0];
rz(-1.8796896) q[0];
rz(-0.17503861) q[1];
sx q[1];
rz(-1.9997528) q[1];
sx q[1];
rz(-1.6040241) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8309098) q[0];
sx q[0];
rz(-1.791782) q[0];
sx q[0];
rz(-1.2587147) q[0];
rz(-pi) q[1];
x q[1];
rz(0.23065718) q[2];
sx q[2];
rz(-2.7129177) q[2];
sx q[2];
rz(0.32282695) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.9201587) q[1];
sx q[1];
rz(-2.4281574) q[1];
sx q[1];
rz(2.9677797) q[1];
rz(2.2213307) q[3];
sx q[3];
rz(-0.55763054) q[3];
sx q[3];
rz(2.695431) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.039915446) q[2];
sx q[2];
rz(-2.1890409) q[2];
sx q[2];
rz(-2.3802479) q[2];
rz(2.2411761) q[3];
sx q[3];
rz(-2.5420928) q[3];
sx q[3];
rz(-3.0925687) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3392357) q[0];
sx q[0];
rz(-0.46827066) q[0];
sx q[0];
rz(-0.21690579) q[0];
rz(-0.63198173) q[1];
sx q[1];
rz(-1.4914373) q[1];
sx q[1];
rz(2.1868618) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2749598) q[0];
sx q[0];
rz(-1.2506335) q[0];
sx q[0];
rz(-0.76036705) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5920792) q[2];
sx q[2];
rz(-1.9351442) q[2];
sx q[2];
rz(3.0669616) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.4090875) q[1];
sx q[1];
rz(-1.115683) q[1];
sx q[1];
rz(-1.6556157) q[1];
x q[2];
rz(-0.93654376) q[3];
sx q[3];
rz(-1.2055956) q[3];
sx q[3];
rz(1.6983502) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.1251936) q[2];
sx q[2];
rz(-1.9026326) q[2];
sx q[2];
rz(0.94474244) q[2];
rz(-2.7567806) q[3];
sx q[3];
rz(-1.1184357) q[3];
sx q[3];
rz(0.95705664) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(-2.5007297) q[0];
sx q[0];
rz(-2.6364115) q[0];
sx q[0];
rz(-1.5873948) q[0];
rz(-2.2568933) q[1];
sx q[1];
rz(-0.90507602) q[1];
sx q[1];
rz(-0.25837635) q[1];
rz(2.5169218) q[2];
sx q[2];
rz(-0.93745898) q[2];
sx q[2];
rz(-2.6425101) q[2];
rz(-0.20529071) q[3];
sx q[3];
rz(-1.1355573) q[3];
sx q[3];
rz(-0.60187403) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
