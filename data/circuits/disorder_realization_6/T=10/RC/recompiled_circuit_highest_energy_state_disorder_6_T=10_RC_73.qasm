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
rz(-2.1498635) q[0];
sx q[0];
rz(-2.5340134) q[0];
sx q[0];
rz(2.1460331) q[0];
rz(0.52752703) q[1];
sx q[1];
rz(-2.1105284) q[1];
sx q[1];
rz(0.48479015) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82370347) q[0];
sx q[0];
rz(-1.4394338) q[0];
sx q[0];
rz(-2.0844368) q[0];
rz(-pi) q[1];
rz(-0.44102168) q[2];
sx q[2];
rz(-0.82953757) q[2];
sx q[2];
rz(-1.5953568) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.165786) q[1];
sx q[1];
rz(-2.1790904) q[1];
sx q[1];
rz(-3.0700141) q[1];
rz(-pi) q[2];
rz(-2.8504478) q[3];
sx q[3];
rz(-1.640794) q[3];
sx q[3];
rz(2.002188) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.506044) q[2];
sx q[2];
rz(-2.1940993) q[2];
sx q[2];
rz(0.49723899) q[2];
rz(-2.5908568) q[3];
sx q[3];
rz(-2.4516055) q[3];
sx q[3];
rz(2.0042888) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15744844) q[0];
sx q[0];
rz(-1.7406311) q[0];
sx q[0];
rz(2.3781811) q[0];
rz(0.58452559) q[1];
sx q[1];
rz(-1.7598563) q[1];
sx q[1];
rz(-1.0113299) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4268734) q[0];
sx q[0];
rz(-1.4397286) q[0];
sx q[0];
rz(0.1161954) q[0];
rz(-pi) q[1];
x q[1];
rz(0.23638921) q[2];
sx q[2];
rz(-1.4589785) q[2];
sx q[2];
rz(0.37416247) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.71537575) q[1];
sx q[1];
rz(-1.9131887) q[1];
sx q[1];
rz(-0.21933783) q[1];
rz(-0.92722757) q[3];
sx q[3];
rz(-0.6112186) q[3];
sx q[3];
rz(-0.54655308) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.6611019) q[2];
sx q[2];
rz(-1.2031518) q[2];
sx q[2];
rz(2.1413546) q[2];
rz(2.5601322) q[3];
sx q[3];
rz(-2.4088819) q[3];
sx q[3];
rz(-1.2015517) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30705273) q[0];
sx q[0];
rz(-0.59891278) q[0];
sx q[0];
rz(2.5839928) q[0];
rz(-2.8343976) q[1];
sx q[1];
rz(-0.59978825) q[1];
sx q[1];
rz(-0.44816005) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4059999) q[0];
sx q[0];
rz(-1.7032529) q[0];
sx q[0];
rz(-2.622582) q[0];
x q[1];
rz(1.2101754) q[2];
sx q[2];
rz(-2.5917987) q[2];
sx q[2];
rz(-0.87042337) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.8338018) q[1];
sx q[1];
rz(-1.543601) q[1];
sx q[1];
rz(-2.2016352) q[1];
rz(0.056164785) q[3];
sx q[3];
rz(-2.7066288) q[3];
sx q[3];
rz(2.5954704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.8223411) q[2];
sx q[2];
rz(-2.4288869) q[2];
sx q[2];
rz(2.5437497) q[2];
rz(-2.6939825) q[3];
sx q[3];
rz(-1.2667344) q[3];
sx q[3];
rz(-3.083057) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7287801) q[0];
sx q[0];
rz(-0.030981177) q[0];
sx q[0];
rz(2.431751) q[0];
rz(2.1624883) q[1];
sx q[1];
rz(-2.282228) q[1];
sx q[1];
rz(-1.3213347) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6431491) q[0];
sx q[0];
rz(-1.6263026) q[0];
sx q[0];
rz(0.70141354) q[0];
x q[1];
rz(0.780402) q[2];
sx q[2];
rz(-1.681904) q[2];
sx q[2];
rz(-1.7826796) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7497341) q[1];
sx q[1];
rz(-0.94288153) q[1];
sx q[1];
rz(0.33291746) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0274326) q[3];
sx q[3];
rz(-1.3001983) q[3];
sx q[3];
rz(2.4222514) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.9909624) q[2];
sx q[2];
rz(-1.7143098) q[2];
sx q[2];
rz(2.4347351) q[2];
rz(-1.7957211) q[3];
sx q[3];
rz(-1.462505) q[3];
sx q[3];
rz(2.4367387) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83955806) q[0];
sx q[0];
rz(-0.78080451) q[0];
sx q[0];
rz(-1.5020405) q[0];
rz(-1.543965) q[1];
sx q[1];
rz(-2.2815956) q[1];
sx q[1];
rz(-1.4942716) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65022138) q[0];
sx q[0];
rz(-1.5691391) q[0];
sx q[0];
rz(-1.0003538) q[0];
rz(0.59539184) q[2];
sx q[2];
rz(-1.5768345) q[2];
sx q[2];
rz(-1.5689107) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.0922904) q[1];
sx q[1];
rz(-2.7100724) q[1];
sx q[1];
rz(0.79705654) q[1];
rz(-pi) q[2];
rz(-3.1190514) q[3];
sx q[3];
rz(-0.95571858) q[3];
sx q[3];
rz(2.817077) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.26022252) q[2];
sx q[2];
rz(-1.2726731) q[2];
sx q[2];
rz(-0.55848813) q[2];
rz(-0.50885606) q[3];
sx q[3];
rz(-1.5285834) q[3];
sx q[3];
rz(-1.8744899) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5452071) q[0];
sx q[0];
rz(-1.243243) q[0];
sx q[0];
rz(2.9811356) q[0];
rz(-2.4682553) q[1];
sx q[1];
rz(-1.9247749) q[1];
sx q[1];
rz(-2.014726) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50235275) q[0];
sx q[0];
rz(-1.9316088) q[0];
sx q[0];
rz(-2.8176184) q[0];
rz(-2.8062988) q[2];
sx q[2];
rz(-0.79293434) q[2];
sx q[2];
rz(-2.0814587) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.7221795) q[1];
sx q[1];
rz(-2.1447138) q[1];
sx q[1];
rz(0.24146059) q[1];
rz(-0.2114919) q[3];
sx q[3];
rz(-1.0387254) q[3];
sx q[3];
rz(-2.718975) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.0141853) q[2];
sx q[2];
rz(-1.5994453) q[2];
sx q[2];
rz(2.3678153) q[2];
rz(0.989178) q[3];
sx q[3];
rz(-2.7343605) q[3];
sx q[3];
rz(1.0729084) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2268552) q[0];
sx q[0];
rz(-1.6586774) q[0];
sx q[0];
rz(-2.3684655) q[0];
rz(1.585656) q[1];
sx q[1];
rz(-1.2890041) q[1];
sx q[1];
rz(1.1770491) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0601412) q[0];
sx q[0];
rz(-0.8340652) q[0];
sx q[0];
rz(0.042515083) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1074064) q[2];
sx q[2];
rz(-2.5167373) q[2];
sx q[2];
rz(1.8996849) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.8428264) q[1];
sx q[1];
rz(-0.93888679) q[1];
sx q[1];
rz(2.7999987) q[1];
rz(0.66399948) q[3];
sx q[3];
rz(-2.917508) q[3];
sx q[3];
rz(-2.1685719) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.3561463) q[2];
sx q[2];
rz(-2.7768504) q[2];
sx q[2];
rz(2.6981603) q[2];
rz(0.38608471) q[3];
sx q[3];
rz(-1.4177136) q[3];
sx q[3];
rz(-2.940787) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3767553) q[0];
sx q[0];
rz(-1.5926462) q[0];
sx q[0];
rz(0.03373294) q[0];
rz(0.69817606) q[1];
sx q[1];
rz(-2.5444784) q[1];
sx q[1];
rz(0.66506344) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3515776) q[0];
sx q[0];
rz(-0.52886334) q[0];
sx q[0];
rz(-2.1109423) q[0];
rz(-pi) q[1];
rz(-2.6680634) q[2];
sx q[2];
rz(-2.2396834) q[2];
sx q[2];
rz(-0.39448002) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.60282502) q[1];
sx q[1];
rz(-2.1884349) q[1];
sx q[1];
rz(-2.5185939) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5499209) q[3];
sx q[3];
rz(-1.890278) q[3];
sx q[3];
rz(-0.11939458) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.4195121) q[2];
sx q[2];
rz(-0.91700143) q[2];
sx q[2];
rz(1.3631932) q[2];
rz(-0.52115399) q[3];
sx q[3];
rz(-1.8173953) q[3];
sx q[3];
rz(0.49191973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1727961) q[0];
sx q[0];
rz(-1.7954614) q[0];
sx q[0];
rz(-0.81934339) q[0];
rz(-0.68483812) q[1];
sx q[1];
rz(-1.2872773) q[1];
sx q[1];
rz(3.0192764) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7193094) q[0];
sx q[0];
rz(-2.3580812) q[0];
sx q[0];
rz(0.95565234) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0367461) q[2];
sx q[2];
rz(-0.66613644) q[2];
sx q[2];
rz(-2.9532847) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.51309359) q[1];
sx q[1];
rz(-2.4446396) q[1];
sx q[1];
rz(-0.51817399) q[1];
rz(-1.1055533) q[3];
sx q[3];
rz(-2.1085582) q[3];
sx q[3];
rz(-0.28013438) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.43716064) q[2];
sx q[2];
rz(-1.8081343) q[2];
sx q[2];
rz(2.2740347) q[2];
rz(1.3982754) q[3];
sx q[3];
rz(-2.7531392) q[3];
sx q[3];
rz(0.5790264) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.687998) q[0];
sx q[0];
rz(-1.2247676) q[0];
sx q[0];
rz(1.0892185) q[0];
rz(3.094063) q[1];
sx q[1];
rz(-2.0952974) q[1];
sx q[1];
rz(-2.4542123) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7318667) q[0];
sx q[0];
rz(-1.0805403) q[0];
sx q[0];
rz(3.114028) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0014761) q[2];
sx q[2];
rz(-1.3146558) q[2];
sx q[2];
rz(1.3821951) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.2694305) q[1];
sx q[1];
rz(-2.3331642) q[1];
sx q[1];
rz(-2.4926659) q[1];
rz(-1.6964558) q[3];
sx q[3];
rz(-2.501125) q[3];
sx q[3];
rz(-2.6844269) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.99531949) q[2];
sx q[2];
rz(-0.61890382) q[2];
sx q[2];
rz(-1.4825561) q[2];
rz(-2.4847374) q[3];
sx q[3];
rz(-1.992179) q[3];
sx q[3];
rz(-0.38203865) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8889846) q[0];
sx q[0];
rz(-1.9883307) q[0];
sx q[0];
rz(-2.0345732) q[0];
rz(-0.57869115) q[1];
sx q[1];
rz(-0.83732579) q[1];
sx q[1];
rz(2.8566828) q[1];
rz(-1.6538229) q[2];
sx q[2];
rz(-1.9258693) q[2];
sx q[2];
rz(-3.0167962) q[2];
rz(-2.6246351) q[3];
sx q[3];
rz(-1.6421972) q[3];
sx q[3];
rz(0.35785892) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
