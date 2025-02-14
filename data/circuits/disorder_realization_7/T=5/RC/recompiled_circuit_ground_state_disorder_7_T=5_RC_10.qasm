OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.6586128) q[0];
sx q[0];
rz(-0.40402544) q[0];
sx q[0];
rz(-2.7513096) q[0];
rz(-0.033493869) q[1];
sx q[1];
rz(-2.610425) q[1];
sx q[1];
rz(0.20294987) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73589486) q[0];
sx q[0];
rz(-0.66758388) q[0];
sx q[0];
rz(2.1632458) q[0];
x q[1];
rz(-3.13596) q[2];
sx q[2];
rz(-1.6597103) q[2];
sx q[2];
rz(2.2064457) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.4642849) q[1];
sx q[1];
rz(-2.5618636) q[1];
sx q[1];
rz(-0.34601684) q[1];
rz(-pi) q[2];
rz(1.6018576) q[3];
sx q[3];
rz(-0.72664562) q[3];
sx q[3];
rz(-2.2769312) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.1927294) q[2];
sx q[2];
rz(-2.296083) q[2];
sx q[2];
rz(1.3884937) q[2];
rz(0.071831547) q[3];
sx q[3];
rz(-2.6303232) q[3];
sx q[3];
rz(0.82740074) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.097505957) q[0];
sx q[0];
rz(-2.9768017) q[0];
sx q[0];
rz(2.6440115) q[0];
rz(1.3961821) q[1];
sx q[1];
rz(-1.0284245) q[1];
sx q[1];
rz(-0.57026774) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3511757) q[0];
sx q[0];
rz(-0.47016682) q[0];
sx q[0];
rz(0.29542653) q[0];
rz(0.79945081) q[2];
sx q[2];
rz(-2.852147) q[2];
sx q[2];
rz(1.67414) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-3.1055821) q[1];
sx q[1];
rz(-2.2911628) q[1];
sx q[1];
rz(1.5237113) q[1];
x q[2];
rz(0.87941283) q[3];
sx q[3];
rz(-2.0447404) q[3];
sx q[3];
rz(2.4963801) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.0294864) q[2];
sx q[2];
rz(-1.6873529) q[2];
sx q[2];
rz(-2.4278329) q[2];
rz(-1.7589689) q[3];
sx q[3];
rz(-2.6191923) q[3];
sx q[3];
rz(-0.4471603) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56938982) q[0];
sx q[0];
rz(-2.1694006) q[0];
sx q[0];
rz(-0.33706459) q[0];
rz(0.49346787) q[1];
sx q[1];
rz(-0.69030535) q[1];
sx q[1];
rz(-0.92672551) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0729869) q[0];
sx q[0];
rz(-1.301501) q[0];
sx q[0];
rz(-0.26871839) q[0];
x q[1];
rz(-0.91442911) q[2];
sx q[2];
rz(-0.85322748) q[2];
sx q[2];
rz(-2.2877483) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.2973916) q[1];
sx q[1];
rz(-1.9866041) q[1];
sx q[1];
rz(-2.9625921) q[1];
x q[2];
rz(-2.6217691) q[3];
sx q[3];
rz(-0.9614203) q[3];
sx q[3];
rz(1.6345616) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.1885163) q[2];
sx q[2];
rz(-0.86543721) q[2];
sx q[2];
rz(0.68149978) q[2];
rz(1.513688) q[3];
sx q[3];
rz(-1.8047921) q[3];
sx q[3];
rz(-2.9708235) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3156768) q[0];
sx q[0];
rz(-0.1033533) q[0];
sx q[0];
rz(2.090825) q[0];
rz(-2.5283165) q[1];
sx q[1];
rz(-0.79137099) q[1];
sx q[1];
rz(1.9176066) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0014718) q[0];
sx q[0];
rz(-1.2419235) q[0];
sx q[0];
rz(-2.4809581) q[0];
rz(0.54398016) q[2];
sx q[2];
rz(-1.0144941) q[2];
sx q[2];
rz(0.73327366) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.5256796) q[1];
sx q[1];
rz(-2.5469994) q[1];
sx q[1];
rz(0.022371624) q[1];
rz(-pi) q[2];
rz(-0.23777407) q[3];
sx q[3];
rz(-2.536155) q[3];
sx q[3];
rz(-1.0501705) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.0497389) q[2];
sx q[2];
rz(-2.3564796) q[2];
sx q[2];
rz(-0.66748691) q[2];
rz(1.5751754) q[3];
sx q[3];
rz(-2.9508041) q[3];
sx q[3];
rz(0.13458399) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62395537) q[0];
sx q[0];
rz(-1.0409545) q[0];
sx q[0];
rz(0.55476302) q[0];
rz(-1.8476716) q[1];
sx q[1];
rz(-0.41176739) q[1];
sx q[1];
rz(-2.256934) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8644476) q[0];
sx q[0];
rz(-2.6106123) q[0];
sx q[0];
rz(-1.5220716) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0366465) q[2];
sx q[2];
rz(-0.92443442) q[2];
sx q[2];
rz(-2.4798648) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.64120871) q[1];
sx q[1];
rz(-0.51189089) q[1];
sx q[1];
rz(-1.1275395) q[1];
x q[2];
rz(-0.44142044) q[3];
sx q[3];
rz(-1.8251287) q[3];
sx q[3];
rz(-2.1770432) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(3.1137997) q[2];
sx q[2];
rz(-2.5774559) q[2];
sx q[2];
rz(-1.0998868) q[2];
rz(-0.22282985) q[3];
sx q[3];
rz(-1.204071) q[3];
sx q[3];
rz(0.13404624) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1564388) q[0];
sx q[0];
rz(-2.8414861) q[0];
sx q[0];
rz(1.1837748) q[0];
rz(1.8249594) q[1];
sx q[1];
rz(-2.1778409) q[1];
sx q[1];
rz(-1.0354985) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4754921) q[0];
sx q[0];
rz(-2.0543092) q[0];
sx q[0];
rz(2.7782281) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.99656099) q[2];
sx q[2];
rz(-0.6960667) q[2];
sx q[2];
rz(1.0889167) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.70724132) q[1];
sx q[1];
rz(-1.7508669) q[1];
sx q[1];
rz(-1.3959914) q[1];
rz(0.76642613) q[3];
sx q[3];
rz(-1.5186678) q[3];
sx q[3];
rz(2.9418569) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.7684795) q[2];
sx q[2];
rz(-0.36102411) q[2];
sx q[2];
rz(0.38144544) q[2];
rz(1.9253731) q[3];
sx q[3];
rz(-0.65665025) q[3];
sx q[3];
rz(2.5372964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
x q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4878047) q[0];
sx q[0];
rz(-2.8233546) q[0];
sx q[0];
rz(0.26746622) q[0];
rz(-1.6194612) q[1];
sx q[1];
rz(-2.5289502) q[1];
sx q[1];
rz(0.47346514) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7077443) q[0];
sx q[0];
rz(-1.1994386) q[0];
sx q[0];
rz(1.9385563) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.47961819) q[2];
sx q[2];
rz(-1.7349744) q[2];
sx q[2];
rz(-2.2043383) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.0021759) q[1];
sx q[1];
rz(-1.3502023) q[1];
sx q[1];
rz(1.6898512) q[1];
rz(0.3977054) q[3];
sx q[3];
rz(-0.4022214) q[3];
sx q[3];
rz(-2.2598852) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.9240616) q[2];
sx q[2];
rz(-1.4950098) q[2];
sx q[2];
rz(-1.7688497) q[2];
rz(-0.30630201) q[3];
sx q[3];
rz(-2.114571) q[3];
sx q[3];
rz(0.33941227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.821625) q[0];
sx q[0];
rz(-0.29682934) q[0];
sx q[0];
rz(-0.4739652) q[0];
rz(0.78556806) q[1];
sx q[1];
rz(-2.5626917) q[1];
sx q[1];
rz(0.89326352) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97800228) q[0];
sx q[0];
rz(-1.8390391) q[0];
sx q[0];
rz(0.017649529) q[0];
x q[1];
rz(-1.1940057) q[2];
sx q[2];
rz(-2.2991141) q[2];
sx q[2];
rz(1.6737936) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.9109237) q[1];
sx q[1];
rz(-1.603447) q[1];
sx q[1];
rz(2.9327675) q[1];
rz(-1.9666709) q[3];
sx q[3];
rz(-1.4787889) q[3];
sx q[3];
rz(2.8407872) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.47024176) q[2];
sx q[2];
rz(-1.525815) q[2];
sx q[2];
rz(2.5435756) q[2];
rz(2.9371069) q[3];
sx q[3];
rz(-3.0308767) q[3];
sx q[3];
rz(2.428875) q[3];
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
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0996899) q[0];
sx q[0];
rz(-2.1429017) q[0];
sx q[0];
rz(0.69910753) q[0];
rz(2.7514669) q[1];
sx q[1];
rz(-2.4603619) q[1];
sx q[1];
rz(-0.9763388) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6481144) q[0];
sx q[0];
rz(-0.29434395) q[0];
sx q[0];
rz(-1.0674547) q[0];
x q[1];
rz(1.2430698) q[2];
sx q[2];
rz(-1.9877745) q[2];
sx q[2];
rz(-2.8627739) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.55667952) q[1];
sx q[1];
rz(-1.0631069) q[1];
sx q[1];
rz(-0.084225144) q[1];
rz(2.5849708) q[3];
sx q[3];
rz(-1.8020013) q[3];
sx q[3];
rz(1.7116261) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.98328996) q[2];
sx q[2];
rz(-0.96902865) q[2];
sx q[2];
rz(2.9950673) q[2];
rz(-2.8807785) q[3];
sx q[3];
rz(-2.0050037) q[3];
sx q[3];
rz(-0.28723106) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89676595) q[0];
sx q[0];
rz(-2.792206) q[0];
sx q[0];
rz(-0.31325999) q[0];
rz(2.3406155) q[1];
sx q[1];
rz(-1.5104834) q[1];
sx q[1];
rz(-0.40447485) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4749194) q[0];
sx q[0];
rz(-1.4316779) q[0];
sx q[0];
rz(0.1864726) q[0];
x q[1];
rz(1.7687665) q[2];
sx q[2];
rz(-2.1083197) q[2];
sx q[2];
rz(-0.24453577) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.7491662) q[1];
sx q[1];
rz(-0.87320864) q[1];
sx q[1];
rz(-2.4568899) q[1];
rz(1.9076212) q[3];
sx q[3];
rz(-1.4416896) q[3];
sx q[3];
rz(2.5548803) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.6910088) q[2];
sx q[2];
rz(-0.49387026) q[2];
sx q[2];
rz(0.79130006) q[2];
rz(2.6386236) q[3];
sx q[3];
rz(-2.0937604) q[3];
sx q[3];
rz(-2.3414229) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7293411) q[0];
sx q[0];
rz(-1.8251735) q[0];
sx q[0];
rz(1.770021) q[0];
rz(-2.5149863) q[1];
sx q[1];
rz(-1.659844) q[1];
sx q[1];
rz(-1.0214092) q[1];
rz(-1.4087497) q[2];
sx q[2];
rz(-1.2037983) q[2];
sx q[2];
rz(-1.498602) q[2];
rz(2.9081591) q[3];
sx q[3];
rz(-1.5830276) q[3];
sx q[3];
rz(2.9434634) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
