OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.72397435) q[0];
sx q[0];
rz(-1.6516049) q[0];
sx q[0];
rz(0.93044257) q[0];
rz(0.62970495) q[1];
sx q[1];
rz(-2.0071109) q[1];
sx q[1];
rz(-1.1073444) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5168092) q[0];
sx q[0];
rz(-0.9073572) q[0];
sx q[0];
rz(-1.9883363) q[0];
rz(-0.8503352) q[2];
sx q[2];
rz(-1.2002266) q[2];
sx q[2];
rz(0.71760273) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.2335637) q[1];
sx q[1];
rz(-1.4587914) q[1];
sx q[1];
rz(-1.2909375) q[1];
rz(-pi) q[2];
rz(-2.2060478) q[3];
sx q[3];
rz(-1.8555292) q[3];
sx q[3];
rz(2.2176544) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.063623108) q[2];
sx q[2];
rz(-2.412553) q[2];
sx q[2];
rz(-1.3280274) q[2];
rz(0.32087457) q[3];
sx q[3];
rz(-0.98595536) q[3];
sx q[3];
rz(3.0096171) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48224738) q[0];
sx q[0];
rz(-0.11238614) q[0];
sx q[0];
rz(-2.2609718) q[0];
rz(-1.847514) q[1];
sx q[1];
rz(-2.7236415) q[1];
sx q[1];
rz(-0.81726384) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0016804455) q[0];
sx q[0];
rz(-1.7117371) q[0];
sx q[0];
rz(-1.9190448) q[0];
x q[1];
rz(1.0802644) q[2];
sx q[2];
rz(-1.3776407) q[2];
sx q[2];
rz(-2.1527388) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.6868999) q[1];
sx q[1];
rz(-0.94140879) q[1];
sx q[1];
rz(1.8886186) q[1];
x q[2];
rz(-1.0975295) q[3];
sx q[3];
rz(-1.739199) q[3];
sx q[3];
rz(-0.67697224) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.2237504) q[2];
sx q[2];
rz(-0.68085256) q[2];
sx q[2];
rz(-2.7775653) q[2];
rz(2.1552127) q[3];
sx q[3];
rz(-1.7247518) q[3];
sx q[3];
rz(-1.4646437) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7746975) q[0];
sx q[0];
rz(-2.3092473) q[0];
sx q[0];
rz(0.96631518) q[0];
rz(0.19293383) q[1];
sx q[1];
rz(-2.0529592) q[1];
sx q[1];
rz(1.4470709) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15196249) q[0];
sx q[0];
rz(-1.785779) q[0];
sx q[0];
rz(1.5887898) q[0];
rz(-0.082712163) q[2];
sx q[2];
rz(-1.2809922) q[2];
sx q[2];
rz(-2.2207584) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.2599064) q[1];
sx q[1];
rz(-0.73598624) q[1];
sx q[1];
rz(-2.5658539) q[1];
rz(-2.0855911) q[3];
sx q[3];
rz(-2.6714973) q[3];
sx q[3];
rz(-1.6197636) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.7029999) q[2];
sx q[2];
rz(-2.6575228) q[2];
sx q[2];
rz(2.036371) q[2];
rz(0.74622074) q[3];
sx q[3];
rz(-1.4893702) q[3];
sx q[3];
rz(2.1658649) q[3];
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
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0063909) q[0];
sx q[0];
rz(-2.6171896) q[0];
sx q[0];
rz(-1.6756469) q[0];
rz(2.8566467) q[1];
sx q[1];
rz(-2.0712712) q[1];
sx q[1];
rz(-2.7526061) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8143512) q[0];
sx q[0];
rz(-1.8530493) q[0];
sx q[0];
rz(-1.3319356) q[0];
x q[1];
rz(-2.94611) q[2];
sx q[2];
rz(-2.0494378) q[2];
sx q[2];
rz(-3.0340956) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.4219141) q[1];
sx q[1];
rz(-0.72243566) q[1];
sx q[1];
rz(-0.6694442) q[1];
rz(1.7763406) q[3];
sx q[3];
rz(-0.83344978) q[3];
sx q[3];
rz(1.8635686) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.4185562) q[2];
sx q[2];
rz(-1.9779466) q[2];
sx q[2];
rz(-2.6848865) q[2];
rz(1.5151954) q[3];
sx q[3];
rz(-0.94907343) q[3];
sx q[3];
rz(-0.28234282) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.221955) q[0];
sx q[0];
rz(-2.0018405) q[0];
sx q[0];
rz(0.36002457) q[0];
rz(-2.4941764) q[1];
sx q[1];
rz(-1.6123687) q[1];
sx q[1];
rz(0.49450758) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3213145) q[0];
sx q[0];
rz(-1.2914133) q[0];
sx q[0];
rz(1.7000291) q[0];
rz(-pi) q[1];
rz(-1.4353384) q[2];
sx q[2];
rz(-1.6098621) q[2];
sx q[2];
rz(-2.1087697) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.20930418) q[1];
sx q[1];
rz(-0.63289019) q[1];
sx q[1];
rz(-3.0401405) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0881808) q[3];
sx q[3];
rz(-2.6969389) q[3];
sx q[3];
rz(-0.85944552) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.71022025) q[2];
sx q[2];
rz(-0.65588313) q[2];
sx q[2];
rz(2.8430856) q[2];
rz(2.8295637) q[3];
sx q[3];
rz(-1.3307064) q[3];
sx q[3];
rz(0.63849866) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1453778) q[0];
sx q[0];
rz(-0.72308102) q[0];
sx q[0];
rz(-0.19590713) q[0];
rz(0.021082489) q[1];
sx q[1];
rz(-1.3985876) q[1];
sx q[1];
rz(-1.9063937) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2081731) q[0];
sx q[0];
rz(-2.6801077) q[0];
sx q[0];
rz(-2.3955406) q[0];
rz(0.70448204) q[2];
sx q[2];
rz(-1.8330169) q[2];
sx q[2];
rz(0.30202497) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.8630353) q[1];
sx q[1];
rz(-2.1457991) q[1];
sx q[1];
rz(-2.0726191) q[1];
rz(-1.487021) q[3];
sx q[3];
rz(-1.3352852) q[3];
sx q[3];
rz(-4/(1*pi)) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6482676) q[2];
sx q[2];
rz(-2.2154634) q[2];
sx q[2];
rz(0.43169272) q[2];
rz(-1.4124983) q[3];
sx q[3];
rz(-2.4192211) q[3];
sx q[3];
rz(3.0055962) q[3];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39111185) q[0];
sx q[0];
rz(-1.9727805) q[0];
sx q[0];
rz(-3.0294763) q[0];
rz(2.926459) q[1];
sx q[1];
rz(-1.5605749) q[1];
sx q[1];
rz(-1.1134061) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8047785) q[0];
sx q[0];
rz(-1.3114197) q[0];
sx q[0];
rz(0.39774261) q[0];
rz(-pi) q[1];
rz(1.3117123) q[2];
sx q[2];
rz(-1.2652745) q[2];
sx q[2];
rz(0.53369001) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.50892936) q[1];
sx q[1];
rz(-2.3718861) q[1];
sx q[1];
rz(2.6973004) q[1];
x q[2];
rz(2.9281499) q[3];
sx q[3];
rz(-1.1323954) q[3];
sx q[3];
rz(-1.160281) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(3.1404184) q[2];
sx q[2];
rz(-0.73264709) q[2];
sx q[2];
rz(0.12602885) q[2];
rz(-2.0942988) q[3];
sx q[3];
rz(-1.8208561) q[3];
sx q[3];
rz(-2.4333911) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
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
rz(0.66013181) q[0];
sx q[0];
rz(-0.75868693) q[0];
sx q[0];
rz(-1.6814167) q[0];
rz(-1.2449645) q[1];
sx q[1];
rz(-1.0943202) q[1];
sx q[1];
rz(-1.2089027) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2395994) q[0];
sx q[0];
rz(-1.308846) q[0];
sx q[0];
rz(-0.082017935) q[0];
x q[1];
rz(-0.65931321) q[2];
sx q[2];
rz(-0.94617832) q[2];
sx q[2];
rz(2.5059932) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.9637451) q[1];
sx q[1];
rz(-0.81804619) q[1];
sx q[1];
rz(-2.6820116) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.21945159) q[3];
sx q[3];
rz(-0.95153522) q[3];
sx q[3];
rz(1.0964364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.37844354) q[2];
sx q[2];
rz(-1.2377137) q[2];
sx q[2];
rz(1.3195999) q[2];
rz(0.59213263) q[3];
sx q[3];
rz(-1.416128) q[3];
sx q[3];
rz(-3.106451) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9086583) q[0];
sx q[0];
rz(-0.52353752) q[0];
sx q[0];
rz(-1.3611025) q[0];
rz(1.9305485) q[1];
sx q[1];
rz(-2.2369604) q[1];
sx q[1];
rz(2.7499054) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2101333) q[0];
sx q[0];
rz(-2.0053232) q[0];
sx q[0];
rz(-1.5643442) q[0];
x q[1];
rz(-2.7462256) q[2];
sx q[2];
rz(-0.48465109) q[2];
sx q[2];
rz(-1.6389099) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.45410941) q[1];
sx q[1];
rz(-0.81067639) q[1];
sx q[1];
rz(2.0337385) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6158995) q[3];
sx q[3];
rz(-1.5075397) q[3];
sx q[3];
rz(-2.4478108) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.52035511) q[2];
sx q[2];
rz(-1.3427799) q[2];
sx q[2];
rz(-1.3367782) q[2];
rz(-2.3796066) q[3];
sx q[3];
rz(-0.31969324) q[3];
sx q[3];
rz(0.80037642) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8005463) q[0];
sx q[0];
rz(-0.30277345) q[0];
sx q[0];
rz(0.57089943) q[0];
rz(1.4292498) q[1];
sx q[1];
rz(-1.0639023) q[1];
sx q[1];
rz(-0.16194078) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8853332) q[0];
sx q[0];
rz(-2.4618039) q[0];
sx q[0];
rz(-1.1023561) q[0];
x q[1];
rz(-3.1324773) q[2];
sx q[2];
rz(-1.7593804) q[2];
sx q[2];
rz(-1.7878704) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.098617741) q[1];
sx q[1];
rz(-0.21511714) q[1];
sx q[1];
rz(0.92623644) q[1];
rz(-pi) q[2];
rz(-1.6041683) q[3];
sx q[3];
rz(-2.0936692) q[3];
sx q[3];
rz(3.0748622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.95742115) q[2];
sx q[2];
rz(-2.0364169) q[2];
sx q[2];
rz(1.7133678) q[2];
rz(1.1994294) q[3];
sx q[3];
rz(-1.0933484) q[3];
sx q[3];
rz(-1.3706346) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4409055) q[0];
sx q[0];
rz(-2.6661243) q[0];
sx q[0];
rz(2.0211924) q[0];
rz(1.7715001) q[1];
sx q[1];
rz(-0.94540989) q[1];
sx q[1];
rz(2.170845) q[1];
rz(1.7182405) q[2];
sx q[2];
rz(-1.4530008) q[2];
sx q[2];
rz(3.0086649) q[2];
rz(0.42967038) q[3];
sx q[3];
rz(-1.3702787) q[3];
sx q[3];
rz(-2.7313781) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
