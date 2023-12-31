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
rz(-2.0071109) q[1];
sx q[1];
rz(-1.1073444) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6247834) q[0];
sx q[0];
rz(-0.9073572) q[0];
sx q[0];
rz(1.1532564) q[0];
x q[1];
rz(2.6644457) q[2];
sx q[2];
rz(-0.90847441) q[2];
sx q[2];
rz(-1.1610111) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.5109374) q[1];
sx q[1];
rz(-1.8488548) q[1];
sx q[1];
rz(-0.11649881) q[1];
x q[2];
rz(-2.2060478) q[3];
sx q[3];
rz(-1.2860635) q[3];
sx q[3];
rz(0.92393827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.063623108) q[2];
sx q[2];
rz(-2.412553) q[2];
sx q[2];
rz(-1.3280274) q[2];
rz(2.8207181) q[3];
sx q[3];
rz(-2.1556373) q[3];
sx q[3];
rz(3.0096171) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6593453) q[0];
sx q[0];
rz(-0.11238614) q[0];
sx q[0];
rz(2.2609718) q[0];
rz(1.847514) q[1];
sx q[1];
rz(-0.41795119) q[1];
sx q[1];
rz(2.3243288) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6200703) q[0];
sx q[0];
rz(-1.9154473) q[0];
sx q[0];
rz(-0.14981139) q[0];
rz(2.9233819) q[2];
sx q[2];
rz(-2.0514224) q[2];
sx q[2];
rz(2.6618119) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.9643147) q[1];
sx q[1];
rz(-0.69523584) q[1];
sx q[1];
rz(0.40538215) q[1];
rz(-pi) q[2];
rz(-1.9278139) q[3];
sx q[3];
rz(-2.6414053) q[3];
sx q[3];
rz(2.564085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.91784224) q[2];
sx q[2];
rz(-2.4607401) q[2];
sx q[2];
rz(2.7775653) q[2];
rz(-2.1552127) q[3];
sx q[3];
rz(-1.7247518) q[3];
sx q[3];
rz(1.4646437) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36689511) q[0];
sx q[0];
rz(-2.3092473) q[0];
sx q[0];
rz(-2.1752775) q[0];
rz(2.9486588) q[1];
sx q[1];
rz(-1.0886334) q[1];
sx q[1];
rz(-1.6945217) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15196249) q[0];
sx q[0];
rz(-1.3558136) q[0];
sx q[0];
rz(-1.5887898) q[0];
rz(-pi) q[1];
rz(-0.082712163) q[2];
sx q[2];
rz(-1.8606004) q[2];
sx q[2];
rz(-0.92083425) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.24070534) q[1];
sx q[1];
rz(-1.9449688) q[1];
sx q[1];
rz(-2.4918873) q[1];
rz(-pi) q[2];
rz(-2.0855911) q[3];
sx q[3];
rz(-0.47009531) q[3];
sx q[3];
rz(-1.521829) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.43859279) q[2];
sx q[2];
rz(-0.48406988) q[2];
sx q[2];
rz(1.1052216) q[2];
rz(2.3953719) q[3];
sx q[3];
rz(-1.6522224) q[3];
sx q[3];
rz(2.1658649) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1352017) q[0];
sx q[0];
rz(-0.52440301) q[0];
sx q[0];
rz(1.6756469) q[0];
rz(-2.8566467) q[1];
sx q[1];
rz(-2.0712712) q[1];
sx q[1];
rz(-0.38898653) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8303191) q[0];
sx q[0];
rz(-1.8000326) q[0];
sx q[0];
rz(-0.29005187) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2126036) q[2];
sx q[2];
rz(-2.6274523) q[2];
sx q[2];
rz(2.6280623) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0491838) q[1];
sx q[1];
rz(-1.0256983) q[1];
sx q[1];
rz(-1.070302) q[1];
rz(-0.74794482) q[3];
sx q[3];
rz(-1.7224632) q[3];
sx q[3];
rz(2.9880854) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.7230364) q[2];
sx q[2];
rz(-1.9779466) q[2];
sx q[2];
rz(0.45670613) q[2];
rz(-1.5151954) q[3];
sx q[3];
rz(-2.1925192) q[3];
sx q[3];
rz(-0.28234282) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
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
rz(-1.6123687) q[1];
sx q[1];
rz(0.49450758) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3213145) q[0];
sx q[0];
rz(-1.2914133) q[0];
sx q[0];
rz(-1.7000291) q[0];
x q[1];
rz(-1.8525271) q[2];
sx q[2];
rz(-3.0006471) q[2];
sx q[2];
rz(0.81705392) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.33489409) q[1];
sx q[1];
rz(-2.1999199) q[1];
sx q[1];
rz(1.6449528) q[1];
rz(1.5962283) q[3];
sx q[3];
rz(-1.126822) q[3];
sx q[3];
rz(-0.80029315) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.71022025) q[2];
sx q[2];
rz(-2.4857095) q[2];
sx q[2];
rz(0.29850706) q[2];
rz(2.8295637) q[3];
sx q[3];
rz(-1.8108862) q[3];
sx q[3];
rz(2.503094) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9962149) q[0];
sx q[0];
rz(-0.72308102) q[0];
sx q[0];
rz(2.9456855) q[0];
rz(-0.021082489) q[1];
sx q[1];
rz(-1.7430051) q[1];
sx q[1];
rz(1.235199) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.93341953) q[0];
sx q[0];
rz(-2.6801077) q[0];
sx q[0];
rz(0.74605201) q[0];
rz(-2.4371106) q[2];
sx q[2];
rz(-1.3085758) q[2];
sx q[2];
rz(2.8395677) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.0683806) q[1];
sx q[1];
rz(-2.3975323) q[1];
sx q[1];
rz(2.5030604) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9052832) q[3];
sx q[3];
rz(-1.4893388) q[3];
sx q[3];
rz(-0.27796516) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.49332508) q[2];
sx q[2];
rz(-0.92612925) q[2];
sx q[2];
rz(0.43169272) q[2];
rz(-1.7290944) q[3];
sx q[3];
rz(-0.72237152) q[3];
sx q[3];
rz(3.0055962) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7504808) q[0];
sx q[0];
rz(-1.9727805) q[0];
sx q[0];
rz(3.0294763) q[0];
rz(-0.21513367) q[1];
sx q[1];
rz(-1.5605749) q[1];
sx q[1];
rz(2.0281866) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3413234) q[0];
sx q[0];
rz(-1.9545262) q[0];
sx q[0];
rz(1.2905489) q[0];
rz(0.31539519) q[2];
sx q[2];
rz(-1.3239667) q[2];
sx q[2];
rz(2.18404) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.6326633) q[1];
sx q[1];
rz(-2.3718861) q[1];
sx q[1];
rz(2.6973004) q[1];
rz(-pi) q[2];
rz(-2.9281499) q[3];
sx q[3];
rz(-2.0091972) q[3];
sx q[3];
rz(-1.160281) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.001174288) q[2];
sx q[2];
rz(-0.73264709) q[2];
sx q[2];
rz(3.0155638) q[2];
rz(-2.0942988) q[3];
sx q[3];
rz(-1.3207366) q[3];
sx q[3];
rz(-0.70820156) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66013181) q[0];
sx q[0];
rz(-0.75868693) q[0];
sx q[0];
rz(1.460176) q[0];
rz(1.8966282) q[1];
sx q[1];
rz(-1.0943202) q[1];
sx q[1];
rz(-1.2089027) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.93223244) q[0];
sx q[0];
rz(-0.27420843) q[0];
sx q[0];
rz(1.8673531) q[0];
rz(0.86645855) q[2];
sx q[2];
rz(-2.2668215) q[2];
sx q[2];
rz(-2.8528086) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.9637451) q[1];
sx q[1];
rz(-0.81804619) q[1];
sx q[1];
rz(-0.45958105) q[1];
rz(1.8672089) q[3];
sx q[3];
rz(-0.65215462) q[3];
sx q[3];
rz(1.46331) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.37844354) q[2];
sx q[2];
rz(-1.903879) q[2];
sx q[2];
rz(-1.8219927) q[2];
rz(-2.54946) q[3];
sx q[3];
rz(-1.416128) q[3];
sx q[3];
rz(0.035141703) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
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
rz(0.2329344) q[0];
sx q[0];
rz(-0.52353752) q[0];
sx q[0];
rz(-1.3611025) q[0];
rz(-1.2110442) q[1];
sx q[1];
rz(-0.90463224) q[1];
sx q[1];
rz(-2.7499054) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4995394) q[0];
sx q[0];
rz(-1.5649438) q[0];
sx q[0];
rz(2.7070579) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6892745) q[2];
sx q[2];
rz(-1.7512133) q[2];
sx q[2];
rz(0.42186055) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.081025) q[1];
sx q[1];
rz(-0.86522663) q[1];
sx q[1];
rz(2.702436) q[1];
x q[2];
rz(-0.12556062) q[3];
sx q[3];
rz(-0.52913044) q[3];
sx q[3];
rz(0.76847968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.6212375) q[2];
sx q[2];
rz(-1.3427799) q[2];
sx q[2];
rz(-1.3367782) q[2];
rz(-0.76198602) q[3];
sx q[3];
rz(-0.31969324) q[3];
sx q[3];
rz(-0.80037642) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8005463) q[0];
sx q[0];
rz(-2.8388192) q[0];
sx q[0];
rz(-2.5706932) q[0];
rz(1.4292498) q[1];
sx q[1];
rz(-1.0639023) q[1];
sx q[1];
rz(-0.16194078) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2019661) q[0];
sx q[0];
rz(-1.2830178) q[0];
sx q[0];
rz(-2.1956325) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3822046) q[2];
sx q[2];
rz(-1.57975) q[2];
sx q[2];
rz(2.9262275) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-3.0429749) q[1];
sx q[1];
rz(-0.21511714) q[1];
sx q[1];
rz(-2.2153562) q[1];
rz(-pi) q[2];
rz(-0.52311388) q[3];
sx q[3];
rz(-1.599708) q[3];
sx q[3];
rz(-1.4873963) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.95742115) q[2];
sx q[2];
rz(-1.1051757) q[2];
sx q[2];
rz(-1.4282248) q[2];
rz(-1.1994294) q[3];
sx q[3];
rz(-1.0933484) q[3];
sx q[3];
rz(1.3706346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7006871) q[0];
sx q[0];
rz(-0.4754684) q[0];
sx q[0];
rz(-1.1204002) q[0];
rz(1.3700925) q[1];
sx q[1];
rz(-2.1961828) q[1];
sx q[1];
rz(-0.97074769) q[1];
rz(-0.8926819) q[2];
sx q[2];
rz(-0.18845367) q[2];
sx q[2];
rz(-1.0343196) q[2];
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
