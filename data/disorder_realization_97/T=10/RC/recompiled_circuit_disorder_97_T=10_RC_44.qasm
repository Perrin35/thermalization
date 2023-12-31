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
rz(-2.2111501) q[0];
rz(0.62970495) q[1];
sx q[1];
rz(4.2760744) q[1];
sx q[1];
rz(8.3174336) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5168092) q[0];
sx q[0];
rz(-2.2342355) q[0];
sx q[0];
rz(1.9883363) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6644457) q[2];
sx q[2];
rz(-2.2331182) q[2];
sx q[2];
rz(1.9805816) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.63065527) q[1];
sx q[1];
rz(-1.2927379) q[1];
sx q[1];
rz(0.11649881) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1125621) q[3];
sx q[3];
rz(-2.4535865) q[3];
sx q[3];
rz(1.0109166) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(3.0779695) q[2];
sx q[2];
rz(-0.72903967) q[2];
sx q[2];
rz(1.8135653) q[2];
rz(-2.8207181) q[3];
sx q[3];
rz(-2.1556373) q[3];
sx q[3];
rz(0.13197556) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48224738) q[0];
sx q[0];
rz(-3.0292065) q[0];
sx q[0];
rz(-2.2609718) q[0];
rz(1.847514) q[1];
sx q[1];
rz(-0.41795119) q[1];
sx q[1];
rz(2.3243288) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0016804455) q[0];
sx q[0];
rz(-1.4298555) q[0];
sx q[0];
rz(1.9190448) q[0];
rz(-pi) q[1];
rz(0.21821071) q[2];
sx q[2];
rz(-1.0901703) q[2];
sx q[2];
rz(-0.47978076) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.9643147) q[1];
sx q[1];
rz(-2.4463568) q[1];
sx q[1];
rz(2.7362105) q[1];
rz(-pi) q[2];
rz(2.9528584) q[3];
sx q[3];
rz(-2.0368324) q[3];
sx q[3];
rz(-2.1621494) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.91784224) q[2];
sx q[2];
rz(-2.4607401) q[2];
sx q[2];
rz(0.36402738) q[2];
rz(-0.98637995) q[3];
sx q[3];
rz(-1.4168408) q[3];
sx q[3];
rz(1.4646437) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7746975) q[0];
sx q[0];
rz(-0.8323454) q[0];
sx q[0];
rz(-2.1752775) q[0];
rz(2.9486588) q[1];
sx q[1];
rz(-1.0886334) q[1];
sx q[1];
rz(1.4470709) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9896302) q[0];
sx q[0];
rz(-1.3558136) q[0];
sx q[0];
rz(-1.5887898) q[0];
rz(0.082712163) q[2];
sx q[2];
rz(-1.2809922) q[2];
sx q[2];
rz(2.2207584) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2599064) q[1];
sx q[1];
rz(-0.73598624) q[1];
sx q[1];
rz(-0.57573872) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0560016) q[3];
sx q[3];
rz(-2.6714973) q[3];
sx q[3];
rz(1.6197636) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.7029999) q[2];
sx q[2];
rz(-2.6575228) q[2];
sx q[2];
rz(1.1052216) q[2];
rz(2.3953719) q[3];
sx q[3];
rz(-1.4893702) q[3];
sx q[3];
rz(0.97572774) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1352017) q[0];
sx q[0];
rz(-0.52440301) q[0];
sx q[0];
rz(1.6756469) q[0];
rz(2.8566467) q[1];
sx q[1];
rz(-2.0712712) q[1];
sx q[1];
rz(-2.7526061) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5327832) q[0];
sx q[0];
rz(-2.7739077) q[0];
sx q[0];
rz(-2.4572548) q[0];
rz(1.2126036) q[2];
sx q[2];
rz(-0.51414031) q[2];
sx q[2];
rz(2.6280623) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.4219141) q[1];
sx q[1];
rz(-0.72243566) q[1];
sx q[1];
rz(2.4721485) q[1];
rz(1.7763406) q[3];
sx q[3];
rz(-0.83344978) q[3];
sx q[3];
rz(-1.2780241) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.4185562) q[2];
sx q[2];
rz(-1.9779466) q[2];
sx q[2];
rz(2.6848865) q[2];
rz(-1.6263973) q[3];
sx q[3];
rz(-2.1925192) q[3];
sx q[3];
rz(0.28234282) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91963768) q[0];
sx q[0];
rz(-1.1397521) q[0];
sx q[0];
rz(2.7815681) q[0];
rz(-2.4941764) q[1];
sx q[1];
rz(-1.5292239) q[1];
sx q[1];
rz(-0.49450758) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8202782) q[0];
sx q[0];
rz(-1.2914133) q[0];
sx q[0];
rz(1.7000291) q[0];
rz(-pi) q[1];
rz(0.039426609) q[2];
sx q[2];
rz(-1.7061503) q[2];
sx q[2];
rz(0.53265041) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.2795909) q[1];
sx q[1];
rz(-1.5108567) q[1];
sx q[1];
rz(0.63043352) q[1];
rz(-pi) q[2];
rz(-0.0534119) q[3];
sx q[3];
rz(-2.6969389) q[3];
sx q[3];
rz(0.85944552) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.4313724) q[2];
sx q[2];
rz(-0.65588313) q[2];
sx q[2];
rz(0.29850706) q[2];
rz(0.31202894) q[3];
sx q[3];
rz(-1.8108862) q[3];
sx q[3];
rz(-2.503094) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9962149) q[0];
sx q[0];
rz(-0.72308102) q[0];
sx q[0];
rz(-2.9456855) q[0];
rz(-0.021082489) q[1];
sx q[1];
rz(-1.7430051) q[1];
sx q[1];
rz(-1.9063937) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7346749) q[0];
sx q[0];
rz(-1.2376681) q[0];
sx q[0];
rz(1.245265) q[0];
rz(-pi) q[1];
x q[1];
rz(0.70448204) q[2];
sx q[2];
rz(-1.8330169) q[2];
sx q[2];
rz(0.30202497) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.27855733) q[1];
sx q[1];
rz(-2.1457991) q[1];
sx q[1];
rz(-1.0689736) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.33552334) q[3];
sx q[3];
rz(-2.8918859) q[3];
sx q[3];
rz(1.6186796) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6482676) q[2];
sx q[2];
rz(-0.92612925) q[2];
sx q[2];
rz(-0.43169272) q[2];
rz(1.4124983) q[3];
sx q[3];
rz(-0.72237152) q[3];
sx q[3];
rz(-0.13599642) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7504808) q[0];
sx q[0];
rz(-1.9727805) q[0];
sx q[0];
rz(0.11211638) q[0];
rz(-0.21513367) q[1];
sx q[1];
rz(-1.5605749) q[1];
sx q[1];
rz(-1.1134061) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3413234) q[0];
sx q[0];
rz(-1.1870664) q[0];
sx q[0];
rz(1.2905489) q[0];
rz(-pi) q[1];
rz(2.4593997) q[2];
sx q[2];
rz(-2.7436514) q[2];
sx q[2];
rz(1.2559909) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.6326633) q[1];
sx q[1];
rz(-2.3718861) q[1];
sx q[1];
rz(2.6973004) q[1];
x q[2];
rz(-1.1464305) q[3];
sx q[3];
rz(-0.48454912) q[3];
sx q[3];
rz(-0.68819203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(3.1404184) q[2];
sx q[2];
rz(-2.4089456) q[2];
sx q[2];
rz(3.0155638) q[2];
rz(1.0472939) q[3];
sx q[3];
rz(-1.3207366) q[3];
sx q[3];
rz(2.4333911) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4814608) q[0];
sx q[0];
rz(-0.75868693) q[0];
sx q[0];
rz(1.460176) q[0];
rz(-1.2449645) q[1];
sx q[1];
rz(-2.0472725) q[1];
sx q[1];
rz(1.2089027) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7891114) q[0];
sx q[0];
rz(-1.6500104) q[0];
sx q[0];
rz(1.8335901) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2751341) q[2];
sx q[2];
rz(-0.87477113) q[2];
sx q[2];
rz(2.8528086) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.066684494) q[1];
sx q[1];
rz(-1.2411331) q[1];
sx q[1];
rz(-2.3782905) q[1];
x q[2];
rz(2.9221411) q[3];
sx q[3];
rz(-2.1900574) q[3];
sx q[3];
rz(-1.0964364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.7631491) q[2];
sx q[2];
rz(-1.2377137) q[2];
sx q[2];
rz(1.3195999) q[2];
rz(2.54946) q[3];
sx q[3];
rz(-1.7254646) q[3];
sx q[3];
rz(-3.106451) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9086583) q[0];
sx q[0];
rz(-2.6180551) q[0];
sx q[0];
rz(1.3611025) q[0];
rz(-1.2110442) q[1];
sx q[1];
rz(-0.90463224) q[1];
sx q[1];
rz(-2.7499054) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4995394) q[0];
sx q[0];
rz(-1.5766489) q[0];
sx q[0];
rz(-0.43453479) q[0];
x q[1];
rz(-1.7708771) q[2];
sx q[2];
rz(-1.1263501) q[2];
sx q[2];
rz(-2.0796298) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.081025) q[1];
sx q[1];
rz(-0.86522663) q[1];
sx q[1];
rz(0.4391567) q[1];
rz(-pi) q[2];
x q[2];
rz(3.016032) q[3];
sx q[3];
rz(-2.6124622) q[3];
sx q[3];
rz(-0.76847968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.52035511) q[2];
sx q[2];
rz(-1.7988127) q[2];
sx q[2];
rz(1.3367782) q[2];
rz(-2.3796066) q[3];
sx q[3];
rz(-0.31969324) q[3];
sx q[3];
rz(0.80037642) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(15/(14*pi)) q[0];
sx q[0];
rz(-2.8388192) q[0];
sx q[0];
rz(2.5706932) q[0];
rz(-1.7123429) q[1];
sx q[1];
rz(-2.0776904) q[1];
sx q[1];
rz(0.16194078) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.93962651) q[0];
sx q[0];
rz(-1.2830178) q[0];
sx q[0];
rz(0.94596011) q[0];
rz(-pi) q[1];
rz(3.1324773) q[2];
sx q[2];
rz(-1.3822123) q[2];
sx q[2];
rz(1.3537223) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.75441277) q[1];
sx q[1];
rz(-1.7422692) q[1];
sx q[1];
rz(-0.13053723) q[1];
rz(0.057823618) q[3];
sx q[3];
rz(-2.6177546) q[3];
sx q[3];
rz(3.0081089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.95742115) q[2];
sx q[2];
rz(-1.1051757) q[2];
sx q[2];
rz(1.4282248) q[2];
rz(-1.9421633) q[3];
sx q[3];
rz(-2.0482443) q[3];
sx q[3];
rz(1.3706346) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4409055) q[0];
sx q[0];
rz(-0.4754684) q[0];
sx q[0];
rz(-1.1204002) q[0];
rz(1.7715001) q[1];
sx q[1];
rz(-0.94540989) q[1];
sx q[1];
rz(2.170845) q[1];
rz(0.11907555) q[2];
sx q[2];
rz(-1.7172114) q[2];
sx q[2];
rz(-1.7211771) q[2];
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
