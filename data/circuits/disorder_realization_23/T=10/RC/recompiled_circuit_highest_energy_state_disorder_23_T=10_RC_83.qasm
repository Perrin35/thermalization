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
rz(0.5956369) q[0];
sx q[0];
rz(-2.73157) q[0];
sx q[0];
rz(-0.79467839) q[0];
rz(-1.0533286) q[1];
sx q[1];
rz(-0.95212189) q[1];
sx q[1];
rz(-1.9917537) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7139706) q[0];
sx q[0];
rz(-0.29610482) q[0];
sx q[0];
rz(1.4480736) q[0];
rz(-pi) q[1];
rz(2.3056411) q[2];
sx q[2];
rz(-1.3516956) q[2];
sx q[2];
rz(2.2676005) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-3.1049624) q[1];
sx q[1];
rz(-1.399002) q[1];
sx q[1];
rz(2.7378332) q[1];
rz(-pi) q[2];
rz(-1.3142228) q[3];
sx q[3];
rz(-2.024641) q[3];
sx q[3];
rz(0.75405771) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.6940234) q[2];
sx q[2];
rz(-2.1489216) q[2];
sx q[2];
rz(0.64219323) q[2];
rz(1.0688952) q[3];
sx q[3];
rz(-1.5690208) q[3];
sx q[3];
rz(1.4890891) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.883413) q[0];
sx q[0];
rz(-0.0029819948) q[0];
sx q[0];
rz(2.6302443) q[0];
rz(1.5072352) q[1];
sx q[1];
rz(-1.8664482) q[1];
sx q[1];
rz(-1.7587597) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-14/(3*pi)) q[0];
sx q[0];
rz(-1.9535747) q[0];
sx q[0];
rz(0.22269999) q[0];
x q[1];
rz(2.7956656) q[2];
sx q[2];
rz(-0.9285766) q[2];
sx q[2];
rz(1.3042892) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9785215) q[1];
sx q[1];
rz(-2.6252665) q[1];
sx q[1];
rz(-2.2341245) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9622127) q[3];
sx q[3];
rz(-0.53032833) q[3];
sx q[3];
rz(-1.8225924) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.6766659) q[2];
sx q[2];
rz(-1.9366465) q[2];
sx q[2];
rz(-0.10291544) q[2];
rz(-2.0788705) q[3];
sx q[3];
rz(-1.1954185) q[3];
sx q[3];
rz(0.11500558) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9677143) q[0];
sx q[0];
rz(-2.102484) q[0];
sx q[0];
rz(-0.1097196) q[0];
rz(-2.8166855) q[1];
sx q[1];
rz(-1.9903245) q[1];
sx q[1];
rz(-0.5330162) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95105458) q[0];
sx q[0];
rz(-1.8673273) q[0];
sx q[0];
rz(-1.7151136) q[0];
x q[1];
rz(-0.052003666) q[2];
sx q[2];
rz(-1.5468165) q[2];
sx q[2];
rz(2.5131186) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.2238206) q[1];
sx q[1];
rz(-0.9448565) q[1];
sx q[1];
rz(-2.2926032) q[1];
rz(-2.0267065) q[3];
sx q[3];
rz(-2.3264936) q[3];
sx q[3];
rz(0.86642735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.816232) q[2];
sx q[2];
rz(-1.8133138) q[2];
sx q[2];
rz(-1.4556966) q[2];
rz(-3.0560737) q[3];
sx q[3];
rz(-2.072008) q[3];
sx q[3];
rz(0.94676179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38086712) q[0];
sx q[0];
rz(-2.4873698) q[0];
sx q[0];
rz(2.4192659) q[0];
rz(-1.5708615) q[1];
sx q[1];
rz(-1.9941565) q[1];
sx q[1];
rz(0.042472366) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6149276) q[0];
sx q[0];
rz(-0.4836429) q[0];
sx q[0];
rz(0.40348224) q[0];
rz(-pi) q[1];
rz(1.4804777) q[2];
sx q[2];
rz(-1.5866163) q[2];
sx q[2];
rz(-2.9243369) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.60085697) q[1];
sx q[1];
rz(-1.3807978) q[1];
sx q[1];
rz(0.79278391) q[1];
rz(-2.9080584) q[3];
sx q[3];
rz(-2.326205) q[3];
sx q[3];
rz(-2.5366207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.2461569) q[2];
sx q[2];
rz(-1.3763206) q[2];
sx q[2];
rz(2.8975272) q[2];
rz(0.79143381) q[3];
sx q[3];
rz(-1.8297198) q[3];
sx q[3];
rz(-0.70146504) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(0.1112082) q[0];
sx q[0];
rz(-0.1143488) q[0];
sx q[0];
rz(-2.9682888) q[0];
rz(2.2068842) q[1];
sx q[1];
rz(-2.2917031) q[1];
sx q[1];
rz(-0.25397837) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7807863) q[0];
sx q[0];
rz(-0.82366952) q[0];
sx q[0];
rz(2.987848) q[0];
rz(-pi) q[1];
x q[1];
rz(0.59400174) q[2];
sx q[2];
rz(-1.8753759) q[2];
sx q[2];
rz(-2.0128743) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.9559814) q[1];
sx q[1];
rz(-2.4479981) q[1];
sx q[1];
rz(-3.1090956) q[1];
rz(-pi) q[2];
rz(-3.0766904) q[3];
sx q[3];
rz(-1.1838473) q[3];
sx q[3];
rz(-2.153819) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.0737334) q[2];
sx q[2];
rz(-2.8974055) q[2];
sx q[2];
rz(-1.7013288) q[2];
rz(0.67617792) q[3];
sx q[3];
rz(-1.22236) q[3];
sx q[3];
rz(-3.0618099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9192231) q[0];
sx q[0];
rz(-1.3105404) q[0];
sx q[0];
rz(2.4134912) q[0];
rz(-1.1657731) q[1];
sx q[1];
rz(-2.246942) q[1];
sx q[1];
rz(0.56070915) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3325719) q[0];
sx q[0];
rz(-2.999358) q[0];
sx q[0];
rz(2.6660772) q[0];
rz(-1.6892904) q[2];
sx q[2];
rz(-2.117422) q[2];
sx q[2];
rz(-2.7989863) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5821895) q[1];
sx q[1];
rz(-1.786937) q[1];
sx q[1];
rz(-1.1690421) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6339763) q[3];
sx q[3];
rz(-1.6681021) q[3];
sx q[3];
rz(-0.41690505) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5460983) q[2];
sx q[2];
rz(-1.4608773) q[2];
sx q[2];
rz(-2.0687436) q[2];
rz(-0.32235518) q[3];
sx q[3];
rz(-2.7287546) q[3];
sx q[3];
rz(2.7505007) q[3];
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
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76674616) q[0];
sx q[0];
rz(-1.640919) q[0];
sx q[0];
rz(0.40819502) q[0];
rz(-2.8152668) q[1];
sx q[1];
rz(-0.34714454) q[1];
sx q[1];
rz(1.4761285) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8360739) q[0];
sx q[0];
rz(-1.3954432) q[0];
sx q[0];
rz(0.039964635) q[0];
rz(1.9397975) q[2];
sx q[2];
rz(-0.59518669) q[2];
sx q[2];
rz(0.068142224) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.0261111) q[1];
sx q[1];
rz(-0.12899765) q[1];
sx q[1];
rz(-2.0088825) q[1];
rz(1.4011151) q[3];
sx q[3];
rz(-1.7998003) q[3];
sx q[3];
rz(-1.218971) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.8866715) q[2];
sx q[2];
rz(-1.6559867) q[2];
sx q[2];
rz(2.2754748) q[2];
rz(-2.4050889) q[3];
sx q[3];
rz(-1.085142) q[3];
sx q[3];
rz(2.2542663) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0199468) q[0];
sx q[0];
rz(-0.045504657) q[0];
sx q[0];
rz(1.2787) q[0];
rz(-2.1389029) q[1];
sx q[1];
rz(-1.8808782) q[1];
sx q[1];
rz(0.444828) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7072413) q[0];
sx q[0];
rz(-2.0929411) q[0];
sx q[0];
rz(2.971071) q[0];
x q[1];
rz(0.72250267) q[2];
sx q[2];
rz(-1.2151657) q[2];
sx q[2];
rz(2.1788545) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.5096566) q[1];
sx q[1];
rz(-2.3705814) q[1];
sx q[1];
rz(1.3994751) q[1];
rz(-2.5505742) q[3];
sx q[3];
rz(-2.6962387) q[3];
sx q[3];
rz(-1.9231918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.7030846) q[2];
sx q[2];
rz(-1.727641) q[2];
sx q[2];
rz(-0.50602305) q[2];
rz(-2.1972726) q[3];
sx q[3];
rz(-0.025886141) q[3];
sx q[3];
rz(0.73616141) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29189062) q[0];
sx q[0];
rz(-2.1493752) q[0];
sx q[0];
rz(-3.131409) q[0];
rz(2.5306375) q[1];
sx q[1];
rz(-0.98038951) q[1];
sx q[1];
rz(-3.0593061) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7381089) q[0];
sx q[0];
rz(-0.62583621) q[0];
sx q[0];
rz(0.66342579) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.23046979) q[2];
sx q[2];
rz(-2.2374212) q[2];
sx q[2];
rz(2.8067592) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.391606) q[1];
sx q[1];
rz(-2.9175903) q[1];
sx q[1];
rz(-2.8620666) q[1];
rz(1.3547586) q[3];
sx q[3];
rz(-2.2240765) q[3];
sx q[3];
rz(-2.9263484) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.33132195) q[2];
sx q[2];
rz(-1.4484582) q[2];
sx q[2];
rz(-2.4324379) q[2];
rz(-1.6592615) q[3];
sx q[3];
rz(-2.5100561) q[3];
sx q[3];
rz(0.46019301) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.822478) q[0];
sx q[0];
rz(-0.8304441) q[0];
sx q[0];
rz(-2.5122232) q[0];
rz(-1.4156226) q[1];
sx q[1];
rz(-1.6547838) q[1];
sx q[1];
rz(-1.9633044) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7348161) q[0];
sx q[0];
rz(-1.9730076) q[0];
sx q[0];
rz(-2.2262385) q[0];
rz(-2.1599102) q[2];
sx q[2];
rz(-0.58393634) q[2];
sx q[2];
rz(-2.6743025) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.32543031) q[1];
sx q[1];
rz(-0.79968444) q[1];
sx q[1];
rz(2.9887448) q[1];
x q[2];
rz(-1.0052106) q[3];
sx q[3];
rz(-2.2095888) q[3];
sx q[3];
rz(-1.8106724) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.87469953) q[2];
sx q[2];
rz(-1.5692254) q[2];
sx q[2];
rz(-2.3229522) q[2];
rz(1.6308174) q[3];
sx q[3];
rz(-1.8418334) q[3];
sx q[3];
rz(0.40217933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64869399) q[0];
sx q[0];
rz(-1.10981) q[0];
sx q[0];
rz(1.4640402) q[0];
rz(-1.2661487) q[1];
sx q[1];
rz(-1.9910973) q[1];
sx q[1];
rz(1.533351) q[1];
rz(1.2513592) q[2];
sx q[2];
rz(-2.4407959) q[2];
sx q[2];
rz(-2.9052225) q[2];
rz(-0.18465445) q[3];
sx q[3];
rz(-1.481537) q[3];
sx q[3];
rz(2.8164345) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
