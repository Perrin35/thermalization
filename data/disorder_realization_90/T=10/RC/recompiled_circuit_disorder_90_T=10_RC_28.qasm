OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.7251627) q[0];
sx q[0];
rz(0.13983146) q[0];
sx q[0];
rz(10.034372) q[0];
rz(0.66863376) q[1];
sx q[1];
rz(-2.2761087) q[1];
sx q[1];
rz(-0.087021526) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.60458175) q[0];
sx q[0];
rz(-0.7365948) q[0];
sx q[0];
rz(-1.0010615) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3743462) q[2];
sx q[2];
rz(-1.478072) q[2];
sx q[2];
rz(0.35284943) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.062292369) q[1];
sx q[1];
rz(-1.0907409) q[1];
sx q[1];
rz(2.9855707) q[1];
x q[2];
rz(-2.7294455) q[3];
sx q[3];
rz(-1.6133568) q[3];
sx q[3];
rz(-1.7736848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(3.0027851) q[2];
sx q[2];
rz(-1.7613208) q[2];
sx q[2];
rz(0.37386093) q[2];
rz(-0.3368245) q[3];
sx q[3];
rz(-1.5954433) q[3];
sx q[3];
rz(-0.22836223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8841298) q[0];
sx q[0];
rz(-0.3733491) q[0];
sx q[0];
rz(1.9467547) q[0];
rz(3.0589814) q[1];
sx q[1];
rz(-1.9742842) q[1];
sx q[1];
rz(3.1412178) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99129358) q[0];
sx q[0];
rz(-0.55643686) q[0];
sx q[0];
rz(-0.87470212) q[0];
rz(-2.2194901) q[2];
sx q[2];
rz(-2.8173335) q[2];
sx q[2];
rz(0.97933724) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.2692711) q[1];
sx q[1];
rz(-1.3355796) q[1];
sx q[1];
rz(2.5260731) q[1];
x q[2];
rz(-0.46044465) q[3];
sx q[3];
rz(-1.5213955) q[3];
sx q[3];
rz(-1.3747665) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3327545) q[2];
sx q[2];
rz(-2.9465582) q[2];
sx q[2];
rz(0.17671281) q[2];
rz(-0.79408944) q[3];
sx q[3];
rz(-0.77787557) q[3];
sx q[3];
rz(0.55364048) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90213838) q[0];
sx q[0];
rz(-2.4376526) q[0];
sx q[0];
rz(-0.40476558) q[0];
rz(1.8602712) q[1];
sx q[1];
rz(-0.57360137) q[1];
sx q[1];
rz(-1.3084897) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.20298) q[0];
sx q[0];
rz(-0.70883195) q[0];
sx q[0];
rz(-0.75074408) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6042125) q[2];
sx q[2];
rz(-1.7747305) q[2];
sx q[2];
rz(-0.98762074) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.1723459) q[1];
sx q[1];
rz(-1.3174651) q[1];
sx q[1];
rz(0.84448703) q[1];
rz(2.951252) q[3];
sx q[3];
rz(-0.44572178) q[3];
sx q[3];
rz(2.1325071) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.9222766) q[2];
sx q[2];
rz(-2.6185161) q[2];
sx q[2];
rz(3.1211839) q[2];
rz(-2.0698047) q[3];
sx q[3];
rz(-2.0139549) q[3];
sx q[3];
rz(-2.4782457) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14169176) q[0];
sx q[0];
rz(-2.2321556) q[0];
sx q[0];
rz(2.2256057) q[0];
rz(0.46332106) q[1];
sx q[1];
rz(-2.0756192) q[1];
sx q[1];
rz(-1.0571009) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6592641) q[0];
sx q[0];
rz(-1.5963012) q[0];
sx q[0];
rz(2.7042424) q[0];
rz(0.82614233) q[2];
sx q[2];
rz(-0.75781265) q[2];
sx q[2];
rz(3.0405424) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6625729) q[1];
sx q[1];
rz(-2.9111007) q[1];
sx q[1];
rz(1.6307994) q[1];
rz(-pi) q[2];
rz(1.3435059) q[3];
sx q[3];
rz(-0.82633457) q[3];
sx q[3];
rz(2.8518761) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.2477734) q[2];
sx q[2];
rz(-0.63988581) q[2];
sx q[2];
rz(-0.34269732) q[2];
rz(1.4438859) q[3];
sx q[3];
rz(-0.87564898) q[3];
sx q[3];
rz(-0.7152344) q[3];
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
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3559568) q[0];
sx q[0];
rz(-2.5781093) q[0];
sx q[0];
rz(-3.0551531) q[0];
rz(1.3899639) q[1];
sx q[1];
rz(-2.2098863) q[1];
sx q[1];
rz(-2.6729029) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2738004) q[0];
sx q[0];
rz(-2.3228354) q[0];
sx q[0];
rz(-1.0206945) q[0];
rz(2.6195171) q[2];
sx q[2];
rz(-0.35208382) q[2];
sx q[2];
rz(-3.0662231) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.3735818) q[1];
sx q[1];
rz(-1.4410767) q[1];
sx q[1];
rz(2.9380161) q[1];
rz(2.3272446) q[3];
sx q[3];
rz(-1.6296248) q[3];
sx q[3];
rz(-1.3853663) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.1420574) q[2];
sx q[2];
rz(-0.88025981) q[2];
sx q[2];
rz(-0.86432499) q[2];
rz(-2.9243829) q[3];
sx q[3];
rz(-0.61621284) q[3];
sx q[3];
rz(2.8488081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42751673) q[0];
sx q[0];
rz(-1.7345411) q[0];
sx q[0];
rz(-2.9445904) q[0];
rz(1.7794094) q[1];
sx q[1];
rz(-1.660659) q[1];
sx q[1];
rz(0.33624712) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67664458) q[0];
sx q[0];
rz(-2.0083545) q[0];
sx q[0];
rz(-2.4462571) q[0];
x q[1];
rz(0.98724987) q[2];
sx q[2];
rz(-2.0835712) q[2];
sx q[2];
rz(0.057904569) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.5832526) q[1];
sx q[1];
rz(-1.5234158) q[1];
sx q[1];
rz(2.5517795) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.398596) q[3];
sx q[3];
rz(-0.88168722) q[3];
sx q[3];
rz(2.8840051) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.6062935) q[2];
sx q[2];
rz(-1.665411) q[2];
sx q[2];
rz(-2.2952648) q[2];
rz(-1.2396631) q[3];
sx q[3];
rz(-1.2318434) q[3];
sx q[3];
rz(-0.0011750778) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3939312) q[0];
sx q[0];
rz(-1.0378391) q[0];
sx q[0];
rz(2.8826662) q[0];
rz(1.7954284) q[1];
sx q[1];
rz(-1.3849473) q[1];
sx q[1];
rz(-1.0940201) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0712229) q[0];
sx q[0];
rz(-2.5734512) q[0];
sx q[0];
rz(2.2413261) q[0];
rz(3.1292186) q[2];
sx q[2];
rz(-1.7933328) q[2];
sx q[2];
rz(2.6518133) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.5309696) q[1];
sx q[1];
rz(-1.5986643) q[1];
sx q[1];
rz(-1.4762957) q[1];
rz(-pi) q[2];
rz(-0.098071531) q[3];
sx q[3];
rz(-1.5489998) q[3];
sx q[3];
rz(1.1636213) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.1743494) q[2];
sx q[2];
rz(-1.7847585) q[2];
sx q[2];
rz(-0.32067498) q[2];
rz(0.43618068) q[3];
sx q[3];
rz(-0.47302055) q[3];
sx q[3];
rz(-0.44803739) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2748579) q[0];
sx q[0];
rz(-0.47645706) q[0];
sx q[0];
rz(-1.0048237) q[0];
rz(0.59016219) q[1];
sx q[1];
rz(-0.90548038) q[1];
sx q[1];
rz(-0.80387962) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.77792149) q[0];
sx q[0];
rz(-1.6102618) q[0];
sx q[0];
rz(-2.5556593) q[0];
rz(1.3857533) q[2];
sx q[2];
rz(-2.1171326) q[2];
sx q[2];
rz(1.1748558) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.56402928) q[1];
sx q[1];
rz(-1.659435) q[1];
sx q[1];
rz(3.1313529) q[1];
rz(-pi) q[2];
x q[2];
rz(0.11532468) q[3];
sx q[3];
rz(-0.66676312) q[3];
sx q[3];
rz(0.61570864) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.33264318) q[2];
sx q[2];
rz(-2.2968764) q[2];
sx q[2];
rz(-2.1311029) q[2];
rz(-2.5381952) q[3];
sx q[3];
rz(-1.5248652) q[3];
sx q[3];
rz(2.5914014) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6036966) q[0];
sx q[0];
rz(-1.2768856) q[0];
sx q[0];
rz(2.7340775) q[0];
rz(0.28911668) q[1];
sx q[1];
rz(-1.1228077) q[1];
sx q[1];
rz(2.3908652) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84839612) q[0];
sx q[0];
rz(-1.8627394) q[0];
sx q[0];
rz(-1.8672724) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4892011) q[2];
sx q[2];
rz(-1.8494693) q[2];
sx q[2];
rz(-2.1747053) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6299906) q[1];
sx q[1];
rz(-1.2641608) q[1];
sx q[1];
rz(2.929045) q[1];
rz(-1.5446072) q[3];
sx q[3];
rz(-1.2695754) q[3];
sx q[3];
rz(-1.6645886) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.0728545) q[2];
sx q[2];
rz(-2.4170503) q[2];
sx q[2];
rz(-1.9753974) q[2];
rz(-1.7769622) q[3];
sx q[3];
rz(-2.3659673) q[3];
sx q[3];
rz(-0.021818074) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5584548) q[0];
sx q[0];
rz(-2.3151509) q[0];
sx q[0];
rz(1.7512084) q[0];
rz(-0.33069262) q[1];
sx q[1];
rz(-0.76534098) q[1];
sx q[1];
rz(-1.6814544) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18171039) q[0];
sx q[0];
rz(-2.070817) q[0];
sx q[0];
rz(2.659003) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0748126) q[2];
sx q[2];
rz(-1.5504312) q[2];
sx q[2];
rz(0.45769826) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.32669386) q[1];
sx q[1];
rz(-1.1732475) q[1];
sx q[1];
rz(-1.3709929) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4246033) q[3];
sx q[3];
rz(-0.90962142) q[3];
sx q[3];
rz(0.17718525) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.34974393) q[2];
sx q[2];
rz(-1.3246374) q[2];
sx q[2];
rz(1.6798518) q[2];
rz(1.1200303) q[3];
sx q[3];
rz(-0.62168613) q[3];
sx q[3];
rz(0.45564836) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6256975) q[0];
sx q[0];
rz(-1.6202171) q[0];
sx q[0];
rz(1.1465999) q[0];
rz(1.3810146) q[1];
sx q[1];
rz(-1.2881423) q[1];
sx q[1];
rz(-1.2013411) q[1];
rz(-2.7881651) q[2];
sx q[2];
rz(-2.3218487) q[2];
sx q[2];
rz(-0.29816366) q[2];
rz(0.67301987) q[3];
sx q[3];
rz(-1.8802079) q[3];
sx q[3];
rz(0.019284266) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];