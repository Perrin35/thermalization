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
rz(-3.0017612) q[0];
sx q[0];
rz(-0.60959417) q[0];
rz(0.66863376) q[1];
sx q[1];
rz(-2.2761087) q[1];
sx q[1];
rz(3.0545711) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.60458175) q[0];
sx q[0];
rz(-2.4049979) q[0];
sx q[0];
rz(2.1405311) q[0];
rz(-0.76724648) q[2];
sx q[2];
rz(-1.6635206) q[2];
sx q[2];
rz(0.35284943) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.39057186) q[1];
sx q[1];
rz(-0.50288548) q[1];
sx q[1];
rz(1.2807756) q[1];
rz(-pi) q[2];
rz(0.41214715) q[3];
sx q[3];
rz(-1.6133568) q[3];
sx q[3];
rz(1.3679078) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.13880754) q[2];
sx q[2];
rz(-1.7613208) q[2];
sx q[2];
rz(2.7677317) q[2];
rz(2.8047681) q[3];
sx q[3];
rz(-1.5954433) q[3];
sx q[3];
rz(-0.22836223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2574629) q[0];
sx q[0];
rz(-2.7682436) q[0];
sx q[0];
rz(-1.194838) q[0];
rz(0.082611235) q[1];
sx q[1];
rz(-1.1673085) q[1];
sx q[1];
rz(-0.00037489051) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9276792) q[0];
sx q[0];
rz(-1.1535026) q[0];
sx q[0];
rz(-0.37950619) q[0];
rz(-0.20034321) q[2];
sx q[2];
rz(-1.8274954) q[2];
sx q[2];
rz(-0.30470195) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.13818571) q[1];
sx q[1];
rz(-0.97461838) q[1];
sx q[1];
rz(-1.8562993) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.11081103) q[3];
sx q[3];
rz(-2.6786945) q[3];
sx q[3];
rz(3.0447931) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.3327545) q[2];
sx q[2];
rz(-0.19503441) q[2];
sx q[2];
rz(2.9648798) q[2];
rz(0.79408944) q[3];
sx q[3];
rz(-0.77787557) q[3];
sx q[3];
rz(2.5879522) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90213838) q[0];
sx q[0];
rz(-2.4376526) q[0];
sx q[0];
rz(-0.40476558) q[0];
rz(-1.8602712) q[1];
sx q[1];
rz(-0.57360137) q[1];
sx q[1];
rz(-1.8331029) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1256516) q[0];
sx q[0];
rz(-1.1106655) q[0];
sx q[0];
rz(-2.5815651) q[0];
x q[1];
rz(1.8070418) q[2];
sx q[2];
rz(-2.0958732) q[2];
sx q[2];
rz(0.46307785) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8150427) q[1];
sx q[1];
rz(-2.3800441) q[1];
sx q[1];
rz(1.1990859) q[1];
rz(-2.951252) q[3];
sx q[3];
rz(-2.6958709) q[3];
sx q[3];
rz(-1.0090855) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.219316) q[2];
sx q[2];
rz(-0.52307659) q[2];
sx q[2];
rz(0.02040872) q[2];
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
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14169176) q[0];
sx q[0];
rz(-2.2321556) q[0];
sx q[0];
rz(-2.2256057) q[0];
rz(2.6782716) q[1];
sx q[1];
rz(-1.0659734) q[1];
sx q[1];
rz(2.0844918) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6592641) q[0];
sx q[0];
rz(-1.5963012) q[0];
sx q[0];
rz(-2.7042424) q[0];
rz(2.5713423) q[2];
sx q[2];
rz(-1.0409365) q[2];
sx q[2];
rz(2.3392764) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.5406487) q[1];
sx q[1];
rz(-1.800866) q[1];
sx q[1];
rz(-3.1275216) q[1];
x q[2];
rz(-0.23987694) q[3];
sx q[3];
rz(-2.3696218) q[3];
sx q[3];
rz(-3.1023657) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.2477734) q[2];
sx q[2];
rz(-0.63988581) q[2];
sx q[2];
rz(2.7988953) q[2];
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
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3559568) q[0];
sx q[0];
rz(-0.56348339) q[0];
sx q[0];
rz(0.086439565) q[0];
rz(1.7516288) q[1];
sx q[1];
rz(-2.2098863) q[1];
sx q[1];
rz(2.6729029) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5422573) q[0];
sx q[0];
rz(-2.242803) q[0];
sx q[0];
rz(0.50962781) q[0];
x q[1];
rz(-2.6195171) q[2];
sx q[2];
rz(-0.35208382) q[2];
sx q[2];
rz(3.0662231) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7680109) q[1];
sx q[1];
rz(-1.7005159) q[1];
sx q[1];
rz(2.9380161) q[1];
rz(-3.0607871) q[3];
sx q[3];
rz(-0.81597933) q[3];
sx q[3];
rz(-2.9007343) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.1420574) q[2];
sx q[2];
rz(-2.2613328) q[2];
sx q[2];
rz(2.2772677) q[2];
rz(0.21720973) q[3];
sx q[3];
rz(-2.5253798) q[3];
sx q[3];
rz(0.29278452) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42751673) q[0];
sx q[0];
rz(-1.4070516) q[0];
sx q[0];
rz(2.9445904) q[0];
rz(-1.7794094) q[1];
sx q[1];
rz(-1.660659) q[1];
sx q[1];
rz(2.8053455) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7173548) q[0];
sx q[0];
rz(-2.3399118) q[0];
sx q[0];
rz(0.63071155) q[0];
x q[1];
rz(-0.7746081) q[2];
sx q[2];
rz(-2.3850072) q[2];
sx q[2];
rz(-2.2677383) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.5583401) q[1];
sx q[1];
rz(-1.6181769) q[1];
sx q[1];
rz(0.58981311) q[1];
x q[2];
rz(1.398596) q[3];
sx q[3];
rz(-2.2599054) q[3];
sx q[3];
rz(2.8840051) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.53529915) q[2];
sx q[2];
rz(-1.665411) q[2];
sx q[2];
rz(2.2952648) q[2];
rz(-1.2396631) q[3];
sx q[3];
rz(-1.9097493) q[3];
sx q[3];
rz(0.0011750778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7476615) q[0];
sx q[0];
rz(-1.0378391) q[0];
sx q[0];
rz(2.8826662) q[0];
rz(1.7954284) q[1];
sx q[1];
rz(-1.3849473) q[1];
sx q[1];
rz(-1.0940201) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0712229) q[0];
sx q[0];
rz(-2.5734512) q[0];
sx q[0];
rz(2.2413261) q[0];
x q[1];
rz(1.5161683) q[2];
sx q[2];
rz(-2.9187181) q[2];
sx q[2];
rz(0.54578997) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.6106231) q[1];
sx q[1];
rz(-1.5429284) q[1];
sx q[1];
rz(1.4762957) q[1];
rz(-2.9225227) q[3];
sx q[3];
rz(-0.1004569) q[3];
sx q[3];
rz(0.62517525) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.9672433) q[2];
sx q[2];
rz(-1.3568342) q[2];
sx q[2];
rz(-2.8209177) q[2];
rz(2.705412) q[3];
sx q[3];
rz(-2.6685721) q[3];
sx q[3];
rz(2.6935553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2748579) q[0];
sx q[0];
rz(-0.47645706) q[0];
sx q[0];
rz(1.0048237) q[0];
rz(-0.59016219) q[1];
sx q[1];
rz(-2.2361123) q[1];
sx q[1];
rz(2.337713) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3748976) q[0];
sx q[0];
rz(-2.156213) q[0];
sx q[0];
rz(-1.5234408) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.55401037) q[2];
sx q[2];
rz(-1.7286574) q[2];
sx q[2];
rz(-2.6487034) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.4623973) q[1];
sx q[1];
rz(-3.052366) q[1];
sx q[1];
rz(-1.4560844) q[1];
rz(-1.4804833) q[3];
sx q[3];
rz(-2.2323425) q[3];
sx q[3];
rz(-2.6722398) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.8089495) q[2];
sx q[2];
rz(-2.2968764) q[2];
sx q[2];
rz(-2.1311029) q[2];
rz(0.60339749) q[3];
sx q[3];
rz(-1.6167275) q[3];
sx q[3];
rz(-2.5914014) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6036966) q[0];
sx q[0];
rz(-1.2768856) q[0];
sx q[0];
rz(0.4075152) q[0];
rz(0.28911668) q[1];
sx q[1];
rz(-2.018785) q[1];
sx q[1];
rz(0.75072748) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84839612) q[0];
sx q[0];
rz(-1.8627394) q[0];
sx q[0];
rz(1.2743203) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4892011) q[2];
sx q[2];
rz(-1.2921234) q[2];
sx q[2];
rz(0.96688731) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.1474485) q[1];
sx q[1];
rz(-1.7732883) q[1];
sx q[1];
rz(1.8840428) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.30131807) q[3];
sx q[3];
rz(-1.5457866) q[3];
sx q[3];
rz(3.0400288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.0728545) q[2];
sx q[2];
rz(-0.72454238) q[2];
sx q[2];
rz(1.1661952) q[2];
rz(-1.7769622) q[3];
sx q[3];
rz(-2.3659673) q[3];
sx q[3];
rz(-0.021818074) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5584548) q[0];
sx q[0];
rz(-2.3151509) q[0];
sx q[0];
rz(-1.7512084) q[0];
rz(-2.8109) q[1];
sx q[1];
rz(-0.76534098) q[1];
sx q[1];
rz(-1.4601382) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6351785) q[0];
sx q[0];
rz(-1.1513452) q[0];
sx q[0];
rz(-2.123453) q[0];
x q[1];
rz(1.6129458) q[2];
sx q[2];
rz(-2.6372006) q[2];
sx q[2];
rz(-1.1500037) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.8148988) q[1];
sx q[1];
rz(-1.1732475) q[1];
sx q[1];
rz(1.3709929) q[1];
x q[2];
rz(-0.7695997) q[3];
sx q[3];
rz(-2.1160612) q[3];
sx q[3];
rz(-1.8850957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.34974393) q[2];
sx q[2];
rz(-1.3246374) q[2];
sx q[2];
rz(1.6798518) q[2];
rz(1.1200303) q[3];
sx q[3];
rz(-2.5199065) q[3];
sx q[3];
rz(2.6859443) q[3];
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
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5158952) q[0];
sx q[0];
rz(-1.6202171) q[0];
sx q[0];
rz(1.1465999) q[0];
rz(1.3810146) q[1];
sx q[1];
rz(-1.2881423) q[1];
sx q[1];
rz(-1.2013411) q[1];
rz(2.3537221) q[2];
sx q[2];
rz(-1.315016) q[2];
sx q[2];
rz(-1.6223326) q[2];
rz(-2.4685728) q[3];
sx q[3];
rz(-1.8802079) q[3];
sx q[3];
rz(0.019284266) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];