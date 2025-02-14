OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.52387828) q[0];
sx q[0];
rz(-2.5472547) q[0];
sx q[0];
rz(-0.30881) q[0];
rz(-2.8438957) q[1];
sx q[1];
rz(-1.2574137) q[1];
sx q[1];
rz(2.5296192) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4966105) q[0];
sx q[0];
rz(-2.5685446) q[0];
sx q[0];
rz(1.7368421) q[0];
rz(-pi) q[1];
rz(2.0711028) q[2];
sx q[2];
rz(-2.8165952) q[2];
sx q[2];
rz(0.18735841) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.6943629) q[1];
sx q[1];
rz(-1.9026533) q[1];
sx q[1];
rz(-2.993078) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4490922) q[3];
sx q[3];
rz(-1.6095265) q[3];
sx q[3];
rz(-2.8722784) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.1674126) q[2];
sx q[2];
rz(-1.1777425) q[2];
sx q[2];
rz(-0.56498945) q[2];
rz(-0.51088339) q[3];
sx q[3];
rz(-2.9027945) q[3];
sx q[3];
rz(-1.5272944) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6086455) q[0];
sx q[0];
rz(-0.2810418) q[0];
sx q[0];
rz(0.99579048) q[0];
rz(2.5536054) q[1];
sx q[1];
rz(-2.7016787) q[1];
sx q[1];
rz(-1.3194552) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.837916) q[0];
sx q[0];
rz(-0.48978031) q[0];
sx q[0];
rz(2.0950926) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2434364) q[2];
sx q[2];
rz(-1.4991949) q[2];
sx q[2];
rz(-2.763235) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.38182575) q[1];
sx q[1];
rz(-0.78227115) q[1];
sx q[1];
rz(1.8657343) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6922582) q[3];
sx q[3];
rz(-2.0503042) q[3];
sx q[3];
rz(3.1288655) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.0413282) q[2];
sx q[2];
rz(-2.1375956) q[2];
sx q[2];
rz(3.1190994) q[2];
rz(-3.0371173) q[3];
sx q[3];
rz(-1.5375117) q[3];
sx q[3];
rz(0.87048602) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81258881) q[0];
sx q[0];
rz(-1.0306232) q[0];
sx q[0];
rz(-1.6382244) q[0];
rz(3.0606048) q[1];
sx q[1];
rz(-2.4826725) q[1];
sx q[1];
rz(1.9714877) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6011172) q[0];
sx q[0];
rz(-0.99760054) q[0];
sx q[0];
rz(-1.3790087) q[0];
rz(-pi) q[1];
rz(-2.4971737) q[2];
sx q[2];
rz(-1.7303409) q[2];
sx q[2];
rz(-1.5089515) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7602424) q[1];
sx q[1];
rz(-2.6406857) q[1];
sx q[1];
rz(1.5111021) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0634138) q[3];
sx q[3];
rz(-1.3510002) q[3];
sx q[3];
rz(-2.6878217) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.66204232) q[2];
sx q[2];
rz(-0.68632555) q[2];
sx q[2];
rz(-3.0983858) q[2];
rz(-0.19733812) q[3];
sx q[3];
rz(-2.1102326) q[3];
sx q[3];
rz(0.432338) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8869121) q[0];
sx q[0];
rz(-2.0557025) q[0];
sx q[0];
rz(-2.8160954) q[0];
rz(-2.2858641) q[1];
sx q[1];
rz(-0.41494644) q[1];
sx q[1];
rz(-1.0861446) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7522858) q[0];
sx q[0];
rz(-2.6323937) q[0];
sx q[0];
rz(0.26637116) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.83502533) q[2];
sx q[2];
rz(-0.97423282) q[2];
sx q[2];
rz(-1.8955829) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.0289183) q[1];
sx q[1];
rz(-0.54593819) q[1];
sx q[1];
rz(1.250099) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7138163) q[3];
sx q[3];
rz(-1.4849097) q[3];
sx q[3];
rz(-0.68616223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.2073652) q[2];
sx q[2];
rz(-0.97565979) q[2];
sx q[2];
rz(-0.18386851) q[2];
rz(-1.3217226) q[3];
sx q[3];
rz(-2.4403641) q[3];
sx q[3];
rz(0.13300657) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4399399) q[0];
sx q[0];
rz(-2.8337605) q[0];
sx q[0];
rz(-2.2046748) q[0];
rz(2.2266455) q[1];
sx q[1];
rz(-0.85275424) q[1];
sx q[1];
rz(1.459704) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2258218) q[0];
sx q[0];
rz(-0.87683319) q[0];
sx q[0];
rz(-3.0208605) q[0];
rz(2.1906846) q[2];
sx q[2];
rz(-1.0703329) q[2];
sx q[2];
rz(-1.9096979) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.062511584) q[1];
sx q[1];
rz(-0.82014232) q[1];
sx q[1];
rz(3.0081577) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7473518) q[3];
sx q[3];
rz(-1.3271204) q[3];
sx q[3];
rz(-1.0765824) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.0978284) q[2];
sx q[2];
rz(-1.8057258) q[2];
sx q[2];
rz(-2.9975927) q[2];
rz(-2.0740267) q[3];
sx q[3];
rz(-0.34393603) q[3];
sx q[3];
rz(2.4501154) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8320134) q[0];
sx q[0];
rz(-2.3342275) q[0];
sx q[0];
rz(-0.76835865) q[0];
rz(0.8575303) q[1];
sx q[1];
rz(-1.0464959) q[1];
sx q[1];
rz(2.5879477) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34617433) q[0];
sx q[0];
rz(-2.6145929) q[0];
sx q[0];
rz(2.0922776) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3411936) q[2];
sx q[2];
rz(-0.78289778) q[2];
sx q[2];
rz(-1.7097434) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.92451292) q[1];
sx q[1];
rz(-1.9929152) q[1];
sx q[1];
rz(-0.98084992) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1646814) q[3];
sx q[3];
rz(-0.94076983) q[3];
sx q[3];
rz(1.5190558) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.386261) q[2];
sx q[2];
rz(-0.43060455) q[2];
sx q[2];
rz(-2.6957896) q[2];
rz(-1.8949932) q[3];
sx q[3];
rz(-1.3695025) q[3];
sx q[3];
rz(2.2915452) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0724532) q[0];
sx q[0];
rz(-0.9599762) q[0];
sx q[0];
rz(0.14007105) q[0];
rz(2.2135997) q[1];
sx q[1];
rz(-1.3368139) q[1];
sx q[1];
rz(1.6915406) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1221913) q[0];
sx q[0];
rz(-0.1941351) q[0];
sx q[0];
rz(0.96786626) q[0];
rz(-2.1004002) q[2];
sx q[2];
rz(-0.43627377) q[2];
sx q[2];
rz(-2.5817007) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.90277744) q[1];
sx q[1];
rz(-1.1401145) q[1];
sx q[1];
rz(2.4611453) q[1];
rz(-pi) q[2];
rz(1.4163903) q[3];
sx q[3];
rz(-1.7503683) q[3];
sx q[3];
rz(1.7056215) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.037584) q[2];
sx q[2];
rz(-0.17825492) q[2];
sx q[2];
rz(3.0577216) q[2];
rz(-2.1508079) q[3];
sx q[3];
rz(-1.1889941) q[3];
sx q[3];
rz(-1.8535463) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24638076) q[0];
sx q[0];
rz(-1.8354494) q[0];
sx q[0];
rz(-0.35686785) q[0];
rz(1.4419979) q[1];
sx q[1];
rz(-2.5394963) q[1];
sx q[1];
rz(3.1325565) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20140245) q[0];
sx q[0];
rz(-1.3497258) q[0];
sx q[0];
rz(1.3371435) q[0];
rz(-pi) q[1];
rz(-2.4414012) q[2];
sx q[2];
rz(-1.8723566) q[2];
sx q[2];
rz(-1.4789326) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.0865514) q[1];
sx q[1];
rz(-2.8541871) q[1];
sx q[1];
rz(1.2432008) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.610454) q[3];
sx q[3];
rz(-1.1278099) q[3];
sx q[3];
rz(-2.7970527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.25702) q[2];
sx q[2];
rz(-1.7947861) q[2];
sx q[2];
rz(2.436077) q[2];
rz(-2.8722615) q[3];
sx q[3];
rz(-2.0651385) q[3];
sx q[3];
rz(2.1721325) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89328289) q[0];
sx q[0];
rz(-0.80535424) q[0];
sx q[0];
rz(-0.70948187) q[0];
rz(2.4995038) q[1];
sx q[1];
rz(-0.44140068) q[1];
sx q[1];
rz(2.716224) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8474397) q[0];
sx q[0];
rz(-2.3340539) q[0];
sx q[0];
rz(2.5320413) q[0];
x q[1];
rz(-2.7284996) q[2];
sx q[2];
rz(-1.0283264) q[2];
sx q[2];
rz(1.422872) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.21103745) q[1];
sx q[1];
rz(-1.1260707) q[1];
sx q[1];
rz(1.6707604) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2669417) q[3];
sx q[3];
rz(-2.1236478) q[3];
sx q[3];
rz(1.0859717) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.8713259) q[2];
sx q[2];
rz(-2.9623803) q[2];
sx q[2];
rz(2.6628009) q[2];
rz(2.791413) q[3];
sx q[3];
rz(-1.1778573) q[3];
sx q[3];
rz(-2.71463) q[3];
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
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66788524) q[0];
sx q[0];
rz(-0.29692867) q[0];
sx q[0];
rz(-2.3042451) q[0];
rz(-0.26652023) q[1];
sx q[1];
rz(-1.8217249) q[1];
sx q[1];
rz(1.7841608) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11649179) q[0];
sx q[0];
rz(-3.0550346) q[0];
sx q[0];
rz(1.5309912) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0783079) q[2];
sx q[2];
rz(-1.6893759) q[2];
sx q[2];
rz(1.103454) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.3597922) q[1];
sx q[1];
rz(-0.6524274) q[1];
sx q[1];
rz(3.1210207) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8164163) q[3];
sx q[3];
rz(-1.1283698) q[3];
sx q[3];
rz(0.57462245) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.9799161) q[2];
sx q[2];
rz(-0.53407532) q[2];
sx q[2];
rz(0.43816379) q[2];
rz(-2.2257889) q[3];
sx q[3];
rz(-1.6180429) q[3];
sx q[3];
rz(2.3384136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
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
rz(2.5220779) q[0];
sx q[0];
rz(-1.2418455) q[0];
sx q[0];
rz(-1.0157304) q[0];
rz(-0.10454128) q[1];
sx q[1];
rz(-1.3019982) q[1];
sx q[1];
rz(-1.7560538) q[1];
rz(1.2404493) q[2];
sx q[2];
rz(-1.6295682) q[2];
sx q[2];
rz(-0.30420797) q[2];
rz(-0.85827479) q[3];
sx q[3];
rz(-1.0714053) q[3];
sx q[3];
rz(-1.3965931) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
