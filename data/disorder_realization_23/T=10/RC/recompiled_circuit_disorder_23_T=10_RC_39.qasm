OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.8795348) q[0];
sx q[0];
rz(-1.4095925) q[0];
sx q[0];
rz(-1.4341266) q[0];
rz(-2.5073476) q[1];
sx q[1];
rz(-0.60159644) q[1];
sx q[1];
rz(-2.7231725) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7941064) q[0];
sx q[0];
rz(-0.54437629) q[0];
sx q[0];
rz(-2.0967336) q[0];
rz(0.11169545) q[2];
sx q[2];
rz(-1.7061491) q[2];
sx q[2];
rz(-2.7370107) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1474255) q[1];
sx q[1];
rz(-0.52417437) q[1];
sx q[1];
rz(2.0736573) q[1];
x q[2];
rz(-1.1316142) q[3];
sx q[3];
rz(-1.3703128) q[3];
sx q[3];
rz(-0.85103121) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.3216386) q[2];
sx q[2];
rz(-1.3691838) q[2];
sx q[2];
rz(0.83797541) q[2];
rz(0.49301246) q[3];
sx q[3];
rz(-2.8686782) q[3];
sx q[3];
rz(-0.078991927) q[3];
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
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3943966) q[0];
sx q[0];
rz(-2.4173739) q[0];
sx q[0];
rz(1.863742) q[0];
rz(-2.9648119) q[1];
sx q[1];
rz(-1.3143833) q[1];
sx q[1];
rz(-0.4321672) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55789253) q[0];
sx q[0];
rz(-1.0520792) q[0];
sx q[0];
rz(1.2516663) q[0];
rz(-pi) q[1];
rz(-1.997666) q[2];
sx q[2];
rz(-1.9383213) q[2];
sx q[2];
rz(1.9232242) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.1951616) q[1];
sx q[1];
rz(-2.6841607) q[1];
sx q[1];
rz(1.3034348) q[1];
x q[2];
rz(1.0057955) q[3];
sx q[3];
rz(-1.5497991) q[3];
sx q[3];
rz(2.3137623) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.5923578) q[2];
sx q[2];
rz(-1.2640415) q[2];
sx q[2];
rz(2.6548927) q[2];
rz(-1.7633847) q[3];
sx q[3];
rz(-1.2599726) q[3];
sx q[3];
rz(-0.53282213) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.10087) q[0];
sx q[0];
rz(-2.3951055) q[0];
sx q[0];
rz(-2.7242463) q[0];
rz(1.6529282) q[1];
sx q[1];
rz(-2.5960943) q[1];
sx q[1];
rz(-0.506385) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34844549) q[0];
sx q[0];
rz(-1.9217102) q[0];
sx q[0];
rz(-2.9532414) q[0];
x q[1];
rz(2.407079) q[2];
sx q[2];
rz(-1.3909512) q[2];
sx q[2];
rz(0.70914662) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.8295146) q[1];
sx q[1];
rz(-1.8784338) q[1];
sx q[1];
rz(2.679146) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9403463) q[3];
sx q[3];
rz(-1.7412211) q[3];
sx q[3];
rz(2.2263262) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.69904077) q[2];
sx q[2];
rz(-2.6802345) q[2];
sx q[2];
rz(0.59147269) q[2];
rz(-2.5555723) q[3];
sx q[3];
rz(-1.932671) q[3];
sx q[3];
rz(1.4311786) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9716924) q[0];
sx q[0];
rz(-0.62830347) q[0];
sx q[0];
rz(-2.3024094) q[0];
rz(-3.1160141) q[1];
sx q[1];
rz(-0.69568101) q[1];
sx q[1];
rz(-1.5930088) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2617333) q[0];
sx q[0];
rz(-1.121959) q[0];
sx q[0];
rz(-3.0405322) q[0];
x q[1];
rz(-0.71009212) q[2];
sx q[2];
rz(-0.93821628) q[2];
sx q[2];
rz(-1.3131504) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.81739391) q[1];
sx q[1];
rz(-2.3014268) q[1];
sx q[1];
rz(-2.6170931) q[1];
rz(0.56460103) q[3];
sx q[3];
rz(-2.6435916) q[3];
sx q[3];
rz(-2.6548487) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.8362391) q[2];
sx q[2];
rz(-0.88399115) q[2];
sx q[2];
rz(3.0419066) q[2];
rz(2.1827407) q[3];
sx q[3];
rz(-1.8226263) q[3];
sx q[3];
rz(1.8166186) q[3];
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
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56617671) q[0];
sx q[0];
rz(-1.7594936) q[0];
sx q[0];
rz(2.8856522) q[0];
rz(0.4610962) q[1];
sx q[1];
rz(-1.0436811) q[1];
sx q[1];
rz(2.3815313) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3292023) q[0];
sx q[0];
rz(-0.15604067) q[0];
sx q[0];
rz(-0.72326707) q[0];
x q[1];
rz(-3.0779482) q[2];
sx q[2];
rz(-1.9846989) q[2];
sx q[2];
rz(-0.49583437) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.8885986) q[1];
sx q[1];
rz(-0.18054403) q[1];
sx q[1];
rz(0.79043364) q[1];
rz(-2.0403967) q[3];
sx q[3];
rz(-2.4393775) q[3];
sx q[3];
rz(-1.5839674) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.5710859) q[2];
sx q[2];
rz(-1.3723624) q[2];
sx q[2];
rz(2.467353) q[2];
rz(-2.9267866) q[3];
sx q[3];
rz(-0.45682296) q[3];
sx q[3];
rz(-0.017344346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53133416) q[0];
sx q[0];
rz(-1.6711618) q[0];
sx q[0];
rz(1.1791139) q[0];
rz(-2.9367661) q[1];
sx q[1];
rz(-2.3463459) q[1];
sx q[1];
rz(-1.0669605) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.897192) q[0];
sx q[0];
rz(-3.1163437) q[0];
sx q[0];
rz(-2.5301754) q[0];
rz(-pi) q[1];
rz(0.42491575) q[2];
sx q[2];
rz(-1.9746466) q[2];
sx q[2];
rz(-1.0645107) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8384335) q[1];
sx q[1];
rz(-1.618297) q[1];
sx q[1];
rz(-2.3179503) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.1596998) q[3];
sx q[3];
rz(-2.3387863) q[3];
sx q[3];
rz(1.476895) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.71010464) q[2];
sx q[2];
rz(-1.8523214) q[2];
sx q[2];
rz(-2.5816494) q[2];
rz(2.4152749) q[3];
sx q[3];
rz(-2.8328219) q[3];
sx q[3];
rz(0.30549756) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2840435) q[0];
sx q[0];
rz(-0.54656583) q[0];
sx q[0];
rz(1.42111) q[0];
rz(0.20206085) q[1];
sx q[1];
rz(-1.7077363) q[1];
sx q[1];
rz(0.85817671) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9439529) q[0];
sx q[0];
rz(-1.4240992) q[0];
sx q[0];
rz(-0.12978817) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0965205) q[2];
sx q[2];
rz(-1.9809234) q[2];
sx q[2];
rz(-2.8571667) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.9550025) q[1];
sx q[1];
rz(-1.5860671) q[1];
sx q[1];
rz(-3.112622) q[1];
rz(-pi) q[2];
rz(1.9167561) q[3];
sx q[3];
rz(-2.1145027) q[3];
sx q[3];
rz(2.791415) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.8537366) q[2];
sx q[2];
rz(-0.47984543) q[2];
sx q[2];
rz(1.8161592) q[2];
rz(-0.89007968) q[3];
sx q[3];
rz(-1.1471014) q[3];
sx q[3];
rz(-0.98852283) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
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
rz(1.6778075) q[0];
sx q[0];
rz(-2.5690434) q[0];
sx q[0];
rz(2.7668787) q[0];
rz(-2.162714) q[1];
sx q[1];
rz(-2.4596877) q[1];
sx q[1];
rz(1.3495061) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9106306) q[0];
sx q[0];
rz(-0.9696784) q[0];
sx q[0];
rz(2.0448951) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4235054) q[2];
sx q[2];
rz(-1.5535083) q[2];
sx q[2];
rz(-1.7503439) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.046557758) q[1];
sx q[1];
rz(-1.3890424) q[1];
sx q[1];
rz(-3.0912116) q[1];
rz(-pi) q[2];
x q[2];
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
rz(0.65138856) q[2];
sx q[2];
rz(-0.75275246) q[2];
sx q[2];
rz(-2.6728969) q[2];
rz(-1.1941341) q[3];
sx q[3];
rz(-1.3041376) q[3];
sx q[3];
rz(-1.3635925) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5114708) q[0];
sx q[0];
rz(-0.90634316) q[0];
sx q[0];
rz(1.8796896) q[0];
rz(-2.966554) q[1];
sx q[1];
rz(-1.9997528) q[1];
sx q[1];
rz(-1.5375686) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14311929) q[0];
sx q[0];
rz(-0.38030085) q[0];
sx q[0];
rz(-2.2024676) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.23065718) q[2];
sx q[2];
rz(-0.42867491) q[2];
sx q[2];
rz(0.32282695) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6919842) q[1];
sx q[1];
rz(-2.2712628) q[1];
sx q[1];
rz(1.4222172) q[1];
x q[2];
rz(0.92026199) q[3];
sx q[3];
rz(-2.5839621) q[3];
sx q[3];
rz(2.695431) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.039915446) q[2];
sx q[2];
rz(-2.1890409) q[2];
sx q[2];
rz(2.3802479) q[2];
rz(-0.90041655) q[3];
sx q[3];
rz(-0.59949985) q[3];
sx q[3];
rz(-0.049023978) q[3];
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
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.802357) q[0];
sx q[0];
rz(-0.46827066) q[0];
sx q[0];
rz(2.9246869) q[0];
rz(-0.63198173) q[1];
sx q[1];
rz(-1.6501553) q[1];
sx q[1];
rz(0.95473081) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.117884) q[0];
sx q[0];
rz(-0.81239359) q[0];
sx q[0];
rz(-2.6931767) q[0];
rz(-pi) q[1];
rz(1.1502613) q[2];
sx q[2];
rz(-1.0610126) q[2];
sx q[2];
rz(-1.860294) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.6001544) q[1];
sx q[1];
rz(-2.6791875) q[1];
sx q[1];
rz(-0.17141436) q[1];
rz(-pi) q[2];
rz(-2.1438164) q[3];
sx q[3];
rz(-2.4224835) q[3];
sx q[3];
rz(-2.5620808) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.1251936) q[2];
sx q[2];
rz(-1.9026326) q[2];
sx q[2];
rz(2.1968502) q[2];
rz(0.38481209) q[3];
sx q[3];
rz(-2.0231569) q[3];
sx q[3];
rz(2.184536) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5007297) q[0];
sx q[0];
rz(-2.6364115) q[0];
sx q[0];
rz(-1.5873948) q[0];
rz(-0.8846994) q[1];
sx q[1];
rz(-2.2365166) q[1];
sx q[1];
rz(2.8832163) q[1];
rz(-0.89818556) q[2];
sx q[2];
rz(-0.85815103) q[2];
sx q[2];
rz(1.382538) q[2];
rz(-1.1273884) q[3];
sx q[3];
rz(-1.3848806) q[3];
sx q[3];
rz(1.0564907) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
