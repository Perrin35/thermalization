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
rz(1.1433831) q[0];
sx q[0];
rz(-1.5902061) q[0];
sx q[0];
rz(-2.8144612) q[0];
rz(-1.3629355) q[1];
sx q[1];
rz(-0.59734806) q[1];
sx q[1];
rz(-2.8391431) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4262213) q[0];
sx q[0];
rz(-1.7914219) q[0];
sx q[0];
rz(2.9721292) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2313868) q[2];
sx q[2];
rz(-1.1288092) q[2];
sx q[2];
rz(-3.0750753) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.55284269) q[1];
sx q[1];
rz(-1.2556013) q[1];
sx q[1];
rz(2.9957459) q[1];
rz(-pi) q[2];
rz(0.20154736) q[3];
sx q[3];
rz(-0.48140654) q[3];
sx q[3];
rz(-2.4498482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.3737619) q[2];
sx q[2];
rz(-0.094002873) q[2];
sx q[2];
rz(2.8685699) q[2];
rz(-2.7772969) q[3];
sx q[3];
rz(-1.7184869) q[3];
sx q[3];
rz(3.1209768) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6881707) q[0];
sx q[0];
rz(-1.0978798) q[0];
sx q[0];
rz(1.706644) q[0];
rz(2.9393328) q[1];
sx q[1];
rz(-0.63019284) q[1];
sx q[1];
rz(0.30466255) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7633491) q[0];
sx q[0];
rz(-1.3225261) q[0];
sx q[0];
rz(2.2909478) q[0];
x q[1];
rz(-2.8424047) q[2];
sx q[2];
rz(-1.9721037) q[2];
sx q[2];
rz(-0.68986675) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.4604291) q[1];
sx q[1];
rz(-0.39084706) q[1];
sx q[1];
rz(2.6397698) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6913139) q[3];
sx q[3];
rz(-1.6778062) q[3];
sx q[3];
rz(-1.4329239) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.27663407) q[2];
sx q[2];
rz(-2.7101597) q[2];
sx q[2];
rz(1.5383833) q[2];
rz(2.3254584) q[3];
sx q[3];
rz(-2.3147801) q[3];
sx q[3];
rz(2.2468467) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64033878) q[0];
sx q[0];
rz(-2.8479939) q[0];
sx q[0];
rz(1.6295992) q[0];
rz(-1.8005796) q[1];
sx q[1];
rz(-0.87313849) q[1];
sx q[1];
rz(0.27892932) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90401201) q[0];
sx q[0];
rz(-1.6991709) q[0];
sx q[0];
rz(-1.6804986) q[0];
rz(-pi) q[1];
x q[1];
rz(0.90143335) q[2];
sx q[2];
rz(-0.7538023) q[2];
sx q[2];
rz(1.19687) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(3.0947436) q[1];
sx q[1];
rz(-1.2654116) q[1];
sx q[1];
rz(2.3161991) q[1];
rz(-pi) q[2];
rz(1.5666714) q[3];
sx q[3];
rz(-0.71625159) q[3];
sx q[3];
rz(0.28820693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.5035847) q[2];
sx q[2];
rz(-1.6562485) q[2];
sx q[2];
rz(0.33835641) q[2];
rz(-1.7342957) q[3];
sx q[3];
rz(-2.0140078) q[3];
sx q[3];
rz(0.97051632) q[3];
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
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19891837) q[0];
sx q[0];
rz(-0.022138683) q[0];
sx q[0];
rz(2.4778147) q[0];
rz(-0.55555073) q[1];
sx q[1];
rz(-1.6352362) q[1];
sx q[1];
rz(1.3551855) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4595918) q[0];
sx q[0];
rz(-0.92456901) q[0];
sx q[0];
rz(-1.610397) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5082325) q[2];
sx q[2];
rz(-2.782397) q[2];
sx q[2];
rz(-2.0421093) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.042204) q[1];
sx q[1];
rz(-0.70800938) q[1];
sx q[1];
rz(0.96570496) q[1];
rz(1.1045214) q[3];
sx q[3];
rz(-1.5111088) q[3];
sx q[3];
rz(0.75327834) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.3163471) q[2];
sx q[2];
rz(-2.1630042) q[2];
sx q[2];
rz(-2.6157731) q[2];
rz(0.60872269) q[3];
sx q[3];
rz(-0.63568297) q[3];
sx q[3];
rz(-0.70613247) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4290001) q[0];
sx q[0];
rz(-1.9569995) q[0];
sx q[0];
rz(-0.96479601) q[0];
rz(-2.7305799) q[1];
sx q[1];
rz(-2.1844468) q[1];
sx q[1];
rz(2.029665) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85025609) q[0];
sx q[0];
rz(-1.5550858) q[0];
sx q[0];
rz(1.5941304) q[0];
rz(-pi) q[1];
rz(-2.0093727) q[2];
sx q[2];
rz(-2.2242332) q[2];
sx q[2];
rz(0.48023047) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.37718038) q[1];
sx q[1];
rz(-2.0736338) q[1];
sx q[1];
rz(-2.2187114) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8741902) q[3];
sx q[3];
rz(-1.7129714) q[3];
sx q[3];
rz(0.09718516) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.3045706) q[2];
sx q[2];
rz(-1.8111753) q[2];
sx q[2];
rz(-2.0988317) q[2];
rz(-1.9581095) q[3];
sx q[3];
rz(-0.8684929) q[3];
sx q[3];
rz(-2.1709501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(0.17827621) q[0];
sx q[0];
rz(-1.9177508) q[0];
sx q[0];
rz(0.29603145) q[0];
rz(-0.0010679642) q[1];
sx q[1];
rz(-1.8031392) q[1];
sx q[1];
rz(-1.0329049) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.293572) q[0];
sx q[0];
rz(-1.5652992) q[0];
sx q[0];
rz(-1.574541) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1994397) q[2];
sx q[2];
rz(-1.6826138) q[2];
sx q[2];
rz(1.5836388) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.34570899) q[1];
sx q[1];
rz(-2.8879037) q[1];
sx q[1];
rz(-2.9182052) q[1];
rz(-pi) q[2];
x q[2];
rz(2.37866) q[3];
sx q[3];
rz(-1.8682533) q[3];
sx q[3];
rz(2.2190588) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.380015) q[2];
sx q[2];
rz(-1.6753316) q[2];
sx q[2];
rz(-1.638691) q[2];
rz(2.5401529) q[3];
sx q[3];
rz(-0.99965874) q[3];
sx q[3];
rz(0.39826605) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0346506) q[0];
sx q[0];
rz(-1.5489738) q[0];
sx q[0];
rz(-0.73367992) q[0];
rz(2.197798) q[1];
sx q[1];
rz(-1.5993886) q[1];
sx q[1];
rz(-0.61818799) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0252903) q[0];
sx q[0];
rz(-2.9553614) q[0];
sx q[0];
rz(1.2982606) q[0];
x q[1];
rz(-1.5603862) q[2];
sx q[2];
rz(-2.4111643) q[2];
sx q[2];
rz(1.324162) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.71431464) q[1];
sx q[1];
rz(-0.66525092) q[1];
sx q[1];
rz(-1.5284431) q[1];
x q[2];
rz(-2.5060081) q[3];
sx q[3];
rz(-1.2762603) q[3];
sx q[3];
rz(2.9377191) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.9435297) q[2];
sx q[2];
rz(-1.1897831) q[2];
sx q[2];
rz(0.61723906) q[2];
rz(-1.9700358) q[3];
sx q[3];
rz(-2.7660683) q[3];
sx q[3];
rz(2.530781) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
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
rz(0.010855762) q[0];
sx q[0];
rz(-2.3541088) q[0];
sx q[0];
rz(2.6686344) q[0];
rz(-1.6821776) q[1];
sx q[1];
rz(-1.4143896) q[1];
sx q[1];
rz(0.032729538) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6115295) q[0];
sx q[0];
rz(-0.034688799) q[0];
sx q[0];
rz(2.7682224) q[0];
rz(2.6647227) q[2];
sx q[2];
rz(-2.3175961) q[2];
sx q[2];
rz(-1.0007953) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(3.1188726) q[1];
sx q[1];
rz(-1.1274844) q[1];
sx q[1];
rz(-2.9037649) q[1];
x q[2];
rz(1.6784689) q[3];
sx q[3];
rz(-1.6009001) q[3];
sx q[3];
rz(-1.6287977) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.1590283) q[2];
sx q[2];
rz(-0.50243598) q[2];
sx q[2];
rz(1.095613) q[2];
rz(-0.76121965) q[3];
sx q[3];
rz(-1.2844362) q[3];
sx q[3];
rz(2.718954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.465076) q[0];
sx q[0];
rz(-3.0618771) q[0];
sx q[0];
rz(-0.7088784) q[0];
rz(-2.830016) q[1];
sx q[1];
rz(-1.0919002) q[1];
sx q[1];
rz(0.30837217) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.209359) q[0];
sx q[0];
rz(-1.4054125) q[0];
sx q[0];
rz(2.7424521) q[0];
rz(-pi) q[1];
rz(-2.7248179) q[2];
sx q[2];
rz(-1.920097) q[2];
sx q[2];
rz(3.0986378) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.31807223) q[1];
sx q[1];
rz(-1.3284577) q[1];
sx q[1];
rz(0.6233539) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3544054) q[3];
sx q[3];
rz(-1.3861674) q[3];
sx q[3];
rz(-2.5291989) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-3.0333905) q[2];
sx q[2];
rz(-1.4806925) q[2];
sx q[2];
rz(-2.9208276) q[2];
rz(-1.0734142) q[3];
sx q[3];
rz(-2.150841) q[3];
sx q[3];
rz(-2.3555135) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49755001) q[0];
sx q[0];
rz(-2.2798517) q[0];
sx q[0];
rz(1.0368689) q[0];
rz(-1.4700302) q[1];
sx q[1];
rz(-2.122888) q[1];
sx q[1];
rz(-1.9482313) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.92591015) q[0];
sx q[0];
rz(-2.0127477) q[0];
sx q[0];
rz(1.6041167) q[0];
rz(0.69865753) q[2];
sx q[2];
rz(-1.9113052) q[2];
sx q[2];
rz(0.27891544) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9876839) q[1];
sx q[1];
rz(-1.7584956) q[1];
sx q[1];
rz(-0.69324268) q[1];
rz(3.0150005) q[3];
sx q[3];
rz(-1.3401573) q[3];
sx q[3];
rz(-2.0960208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.4736453) q[2];
sx q[2];
rz(-1.8569943) q[2];
sx q[2];
rz(3.0038707) q[2];
rz(-1.5406476) q[3];
sx q[3];
rz(-2.9100304) q[3];
sx q[3];
rz(0.93129492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4731428) q[0];
sx q[0];
rz(-1.6669597) q[0];
sx q[0];
rz(-1.054635) q[0];
rz(-0.4857546) q[1];
sx q[1];
rz(-1.2243441) q[1];
sx q[1];
rz(-1.4996554) q[1];
rz(-0.29781801) q[2];
sx q[2];
rz(-0.70651307) q[2];
sx q[2];
rz(2.1950864) q[2];
rz(-1.0578591) q[3];
sx q[3];
rz(-2.2649962) q[3];
sx q[3];
rz(-1.0640127) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
