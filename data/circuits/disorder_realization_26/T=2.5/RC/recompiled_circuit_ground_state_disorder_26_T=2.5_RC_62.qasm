OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.5717918) q[0];
sx q[0];
rz(-0.90919149) q[0];
sx q[0];
rz(-1.6348913) q[0];
rz(1.2070967) q[1];
sx q[1];
rz(6.7758898) q[1];
sx q[1];
rz(13.125782) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6772422) q[0];
sx q[0];
rz(-1.4046245) q[0];
sx q[0];
rz(0.31810036) q[0];
rz(1.1790397) q[2];
sx q[2];
rz(-1.9843352) q[2];
sx q[2];
rz(2.3512083) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.8313347) q[1];
sx q[1];
rz(-1.2911011) q[1];
sx q[1];
rz(-3.1344942) q[1];
rz(-pi) q[2];
rz(2.2484965) q[3];
sx q[3];
rz(-2.428125) q[3];
sx q[3];
rz(1.5763855) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.59014615) q[2];
sx q[2];
rz(-0.85942736) q[2];
sx q[2];
rz(2.0476511) q[2];
rz(-1.236773) q[3];
sx q[3];
rz(-1.1636795) q[3];
sx q[3];
rz(2.8804603) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7673489) q[0];
sx q[0];
rz(-1.2232895) q[0];
sx q[0];
rz(0.71794024) q[0];
rz(-1.9473437) q[1];
sx q[1];
rz(-0.4078882) q[1];
sx q[1];
rz(1.6494707) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30129967) q[0];
sx q[0];
rz(-1.7027149) q[0];
sx q[0];
rz(1.6961845) q[0];
rz(-pi) q[1];
rz(0.43955438) q[2];
sx q[2];
rz(-1.2258523) q[2];
sx q[2];
rz(-2.3203691) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.34864214) q[1];
sx q[1];
rz(-1.1224282) q[1];
sx q[1];
rz(-2.072529) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3371643) q[3];
sx q[3];
rz(-2.3827887) q[3];
sx q[3];
rz(-0.6060588) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.0481999) q[2];
sx q[2];
rz(-0.77646774) q[2];
sx q[2];
rz(-0.36273599) q[2];
rz(-0.87099606) q[3];
sx q[3];
rz(-1.3716776) q[3];
sx q[3];
rz(-1.3236275) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0469623) q[0];
sx q[0];
rz(-1.8383263) q[0];
sx q[0];
rz(0.56075019) q[0];
rz(0.97460711) q[1];
sx q[1];
rz(-2.2798996) q[1];
sx q[1];
rz(-2.9240756) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1851674) q[0];
sx q[0];
rz(-2.5568058) q[0];
sx q[0];
rz(3.1068012) q[0];
rz(-pi) q[1];
rz(1.0571805) q[2];
sx q[2];
rz(-1.4808146) q[2];
sx q[2];
rz(-1.0124026) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.2524061) q[1];
sx q[1];
rz(-1.4562618) q[1];
sx q[1];
rz(-1.3088398) q[1];
rz(-0.09017011) q[3];
sx q[3];
rz(-0.61093447) q[3];
sx q[3];
rz(2.4028525) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.64323032) q[2];
sx q[2];
rz(-2.3537894) q[2];
sx q[2];
rz(-2.1882449) q[2];
rz(2.0638454) q[3];
sx q[3];
rz(-1.6809623) q[3];
sx q[3];
rz(0.47422153) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0039907) q[0];
sx q[0];
rz(-0.1739665) q[0];
sx q[0];
rz(0.91651383) q[0];
rz(0.42463955) q[1];
sx q[1];
rz(-1.3820796) q[1];
sx q[1];
rz(-1.128208) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71867585) q[0];
sx q[0];
rz(-2.2683218) q[0];
sx q[0];
rz(-0.6573702) q[0];
rz(-1.7776971) q[2];
sx q[2];
rz(-0.82447663) q[2];
sx q[2];
rz(-1.4109703) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.27552106) q[1];
sx q[1];
rz(-1.0539319) q[1];
sx q[1];
rz(1.2244774) q[1];
x q[2];
rz(1.8966394) q[3];
sx q[3];
rz(-0.62281194) q[3];
sx q[3];
rz(-2.4145221) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8757223) q[2];
sx q[2];
rz(-1.46572) q[2];
sx q[2];
rz(-1.9423368) q[2];
rz(-1.1572329) q[3];
sx q[3];
rz(-0.80409378) q[3];
sx q[3];
rz(0.51586241) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1542094) q[0];
sx q[0];
rz(-3.0976384) q[0];
sx q[0];
rz(-2.2162345) q[0];
rz(-2.5788653) q[1];
sx q[1];
rz(-1.0147164) q[1];
sx q[1];
rz(2.7664807) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59835382) q[0];
sx q[0];
rz(-2.1133917) q[0];
sx q[0];
rz(1.6981237) q[0];
rz(-pi) q[1];
rz(0.053370287) q[2];
sx q[2];
rz(-1.1359906) q[2];
sx q[2];
rz(-2.82719) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.8567849) q[1];
sx q[1];
rz(-1.4018267) q[1];
sx q[1];
rz(1.595975) q[1];
rz(1.1697606) q[3];
sx q[3];
rz(-1.1479064) q[3];
sx q[3];
rz(2.324375) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.90317) q[2];
sx q[2];
rz(-0.91707245) q[2];
sx q[2];
rz(2.6440716) q[2];
rz(0.43429747) q[3];
sx q[3];
rz(-1.6853761) q[3];
sx q[3];
rz(-2.1527877) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76181805) q[0];
sx q[0];
rz(-0.86123818) q[0];
sx q[0];
rz(-2.8616943) q[0];
rz(2.1474536) q[1];
sx q[1];
rz(-0.86052624) q[1];
sx q[1];
rz(0.57847374) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29541649) q[0];
sx q[0];
rz(-0.65447411) q[0];
sx q[0];
rz(-0.64864463) q[0];
rz(-pi) q[1];
rz(-0.45409338) q[2];
sx q[2];
rz(-0.75579772) q[2];
sx q[2];
rz(2.2735655) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.0117409) q[1];
sx q[1];
rz(-1.3756344) q[1];
sx q[1];
rz(0.73634781) q[1];
rz(-0.26831823) q[3];
sx q[3];
rz(-1.6683516) q[3];
sx q[3];
rz(-0.25948157) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.0591187) q[2];
sx q[2];
rz(-2.5422577) q[2];
sx q[2];
rz(2.2590051) q[2];
rz(-2.5146218) q[3];
sx q[3];
rz(-2.4866703) q[3];
sx q[3];
rz(2.2466808) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51666981) q[0];
sx q[0];
rz(-2.1997917) q[0];
sx q[0];
rz(3.140977) q[0];
rz(-0.90283886) q[1];
sx q[1];
rz(-0.87264624) q[1];
sx q[1];
rz(0.045086233) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0233399) q[0];
sx q[0];
rz(-1.0824507) q[0];
sx q[0];
rz(1.6581384) q[0];
rz(2.9125764) q[2];
sx q[2];
rz(-0.75521246) q[2];
sx q[2];
rz(1.7984185) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.871094) q[1];
sx q[1];
rz(-1.2648858) q[1];
sx q[1];
rz(-2.3853777) q[1];
x q[2];
rz(0.99211971) q[3];
sx q[3];
rz(-1.0401772) q[3];
sx q[3];
rz(1.9506426) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.8467329) q[2];
sx q[2];
rz(-2.5856057) q[2];
sx q[2];
rz(0.9168469) q[2];
rz(0.89546853) q[3];
sx q[3];
rz(-1.6296891) q[3];
sx q[3];
rz(1.2112613) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.004892) q[0];
sx q[0];
rz(-3.0558375) q[0];
sx q[0];
rz(2.6766747) q[0];
rz(1.9526941) q[1];
sx q[1];
rz(-1.7949972) q[1];
sx q[1];
rz(-2.7704923) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44921267) q[0];
sx q[0];
rz(-1.1638196) q[0];
sx q[0];
rz(-2.0203585) q[0];
rz(-pi) q[1];
rz(-3.1041652) q[2];
sx q[2];
rz(-0.82035318) q[2];
sx q[2];
rz(-1.6714753) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.5828723) q[1];
sx q[1];
rz(-2.0077188) q[1];
sx q[1];
rz(-2.3135499) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5072672) q[3];
sx q[3];
rz(-2.6765569) q[3];
sx q[3];
rz(0.1022235) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.4798639) q[2];
sx q[2];
rz(-0.78833818) q[2];
sx q[2];
rz(2.7313477) q[2];
rz(1.6346301) q[3];
sx q[3];
rz(-2.7362636) q[3];
sx q[3];
rz(0.043225616) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.048454408) q[0];
sx q[0];
rz(-0.52670902) q[0];
sx q[0];
rz(2.9658537) q[0];
rz(-2.3161092) q[1];
sx q[1];
rz(-1.8800507) q[1];
sx q[1];
rz(0.82297355) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0153962) q[0];
sx q[0];
rz(-2.1103729) q[0];
sx q[0];
rz(0.71643512) q[0];
rz(0.94613593) q[2];
sx q[2];
rz(-0.51057928) q[2];
sx q[2];
rz(-2.7522786) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.8205679) q[1];
sx q[1];
rz(-1.9181644) q[1];
sx q[1];
rz(0.72183164) q[1];
rz(1.4863344) q[3];
sx q[3];
rz(-2.7296187) q[3];
sx q[3];
rz(2.1503445) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.2814111) q[2];
sx q[2];
rz(-1.1450291) q[2];
sx q[2];
rz(2.6004876) q[2];
rz(-0.42630729) q[3];
sx q[3];
rz(-2.6796894) q[3];
sx q[3];
rz(0.84686744) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4538039) q[0];
sx q[0];
rz(-2.332088) q[0];
sx q[0];
rz(2.8059106) q[0];
rz(3.1188534) q[1];
sx q[1];
rz(-0.72240654) q[1];
sx q[1];
rz(-0.67363277) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1789933) q[0];
sx q[0];
rz(-0.99267753) q[0];
sx q[0];
rz(-1.9890673) q[0];
x q[1];
rz(-2.6725476) q[2];
sx q[2];
rz(-1.6613591) q[2];
sx q[2];
rz(0.064571206) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.4810197) q[1];
sx q[1];
rz(-2.5169454) q[1];
sx q[1];
rz(-1.2195132) q[1];
rz(-0.38506759) q[3];
sx q[3];
rz(-0.35849471) q[3];
sx q[3];
rz(2.8733503) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.444904) q[2];
sx q[2];
rz(-0.99385571) q[2];
sx q[2];
rz(-2.135684) q[2];
rz(2.3575947) q[3];
sx q[3];
rz(-0.32117143) q[3];
sx q[3];
rz(1.706749) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1226596) q[0];
sx q[0];
rz(-1.4228595) q[0];
sx q[0];
rz(2.0056437) q[0];
rz(2.7574273) q[1];
sx q[1];
rz(-1.9728248) q[1];
sx q[1];
rz(-1.493175) q[1];
rz(-0.21239077) q[2];
sx q[2];
rz(-1.3708541) q[2];
sx q[2];
rz(-0.41920991) q[2];
rz(-0.85687153) q[3];
sx q[3];
rz(-1.2752493) q[3];
sx q[3];
rz(-2.8528438) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
