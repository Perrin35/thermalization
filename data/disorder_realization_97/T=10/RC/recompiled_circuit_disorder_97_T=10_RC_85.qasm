OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.4176183) q[0];
sx q[0];
rz(-1.4899878) q[0];
sx q[0];
rz(2.2111501) q[0];
rz(0.62970495) q[1];
sx q[1];
rz(-2.0071109) q[1];
sx q[1];
rz(-1.1073444) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0005172) q[0];
sx q[0];
rz(-0.76671769) q[0];
sx q[0];
rz(-0.47857743) q[0];
rz(-pi) q[1];
rz(-0.8503352) q[2];
sx q[2];
rz(-1.9413661) q[2];
sx q[2];
rz(-0.71760273) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.90802898) q[1];
sx q[1];
rz(-1.4587914) q[1];
sx q[1];
rz(1.2909375) q[1];
rz(-pi) q[2];
rz(-1.1125621) q[3];
sx q[3];
rz(-0.68800612) q[3];
sx q[3];
rz(-2.1306761) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.063623108) q[2];
sx q[2];
rz(-0.72903967) q[2];
sx q[2];
rz(1.8135653) q[2];
rz(0.32087457) q[3];
sx q[3];
rz(-2.1556373) q[3];
sx q[3];
rz(-3.0096171) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48224738) q[0];
sx q[0];
rz(-0.11238614) q[0];
sx q[0];
rz(2.2609718) q[0];
rz(-1.847514) q[1];
sx q[1];
rz(-0.41795119) q[1];
sx q[1];
rz(-2.3243288) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0016804455) q[0];
sx q[0];
rz(-1.4298555) q[0];
sx q[0];
rz(1.9190448) q[0];
rz(-pi) q[1];
rz(-2.9233819) q[2];
sx q[2];
rz(-2.0514224) q[2];
sx q[2];
rz(-2.6618119) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.0664132) q[1];
sx q[1];
rz(-1.8261837) q[1];
sx q[1];
rz(2.4875719) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9278139) q[3];
sx q[3];
rz(-2.6414053) q[3];
sx q[3];
rz(-2.564085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.2237504) q[2];
sx q[2];
rz(-0.68085256) q[2];
sx q[2];
rz(-0.36402738) q[2];
rz(2.1552127) q[3];
sx q[3];
rz(-1.4168408) q[3];
sx q[3];
rz(-1.6769489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36689511) q[0];
sx q[0];
rz(-0.8323454) q[0];
sx q[0];
rz(2.1752775) q[0];
rz(0.19293383) q[1];
sx q[1];
rz(-2.0529592) q[1];
sx q[1];
rz(1.4470709) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15196249) q[0];
sx q[0];
rz(-1.785779) q[0];
sx q[0];
rz(-1.5887898) q[0];
rz(-0.082712163) q[2];
sx q[2];
rz(-1.2809922) q[2];
sx q[2];
rz(0.92083425) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.88168624) q[1];
sx q[1];
rz(-0.73598624) q[1];
sx q[1];
rz(0.57573872) q[1];
x q[2];
rz(-0.2451285) q[3];
sx q[3];
rz(-1.9760625) q[3];
sx q[3];
rz(2.1851636) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7029999) q[2];
sx q[2];
rz(-2.6575228) q[2];
sx q[2];
rz(2.036371) q[2];
rz(-0.74622074) q[3];
sx q[3];
rz(-1.6522224) q[3];
sx q[3];
rz(-0.97572774) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1352017) q[0];
sx q[0];
rz(-0.52440301) q[0];
sx q[0];
rz(-1.6756469) q[0];
rz(2.8566467) q[1];
sx q[1];
rz(-1.0703215) q[1];
sx q[1];
rz(-0.38898653) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31127351) q[0];
sx q[0];
rz(-1.3415601) q[0];
sx q[0];
rz(2.8515408) q[0];
x q[1];
rz(2.94611) q[2];
sx q[2];
rz(-2.0494378) q[2];
sx q[2];
rz(3.0340956) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.3868689) q[1];
sx q[1];
rz(-1.9935973) q[1];
sx q[1];
rz(0.60476426) q[1];
rz(-pi) q[2];
rz(2.9205434) q[3];
sx q[3];
rz(-2.3813558) q[3];
sx q[3];
rz(-1.5628712) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.7230364) q[2];
sx q[2];
rz(-1.1636461) q[2];
sx q[2];
rz(2.6848865) q[2];
rz(-1.5151954) q[3];
sx q[3];
rz(-0.94907343) q[3];
sx q[3];
rz(0.28234282) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91963768) q[0];
sx q[0];
rz(-2.0018405) q[0];
sx q[0];
rz(0.36002457) q[0];
rz(2.4941764) q[1];
sx q[1];
rz(-1.6123687) q[1];
sx q[1];
rz(-0.49450758) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3798824) q[0];
sx q[0];
rz(-2.8344791) q[0];
sx q[0];
rz(2.7193927) q[0];
rz(-pi) q[1];
x q[1];
rz(0.039426609) q[2];
sx q[2];
rz(-1.7061503) q[2];
sx q[2];
rz(-2.6089422) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.20930418) q[1];
sx q[1];
rz(-0.63289019) q[1];
sx q[1];
rz(3.0401405) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5453643) q[3];
sx q[3];
rz(-2.0147707) q[3];
sx q[3];
rz(-2.3412995) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.4313724) q[2];
sx q[2];
rz(-0.65588313) q[2];
sx q[2];
rz(2.8430856) q[2];
rz(0.31202894) q[3];
sx q[3];
rz(-1.3307064) q[3];
sx q[3];
rz(2.503094) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1453778) q[0];
sx q[0];
rz(-0.72308102) q[0];
sx q[0];
rz(0.19590713) q[0];
rz(-3.1205102) q[1];
sx q[1];
rz(-1.7430051) q[1];
sx q[1];
rz(-1.235199) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7346749) q[0];
sx q[0];
rz(-1.2376681) q[0];
sx q[0];
rz(1.8963277) q[0];
x q[1];
rz(0.70448204) q[2];
sx q[2];
rz(-1.3085758) q[2];
sx q[2];
rz(-0.30202497) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.8630353) q[1];
sx q[1];
rz(-0.99579358) q[1];
sx q[1];
rz(1.0689736) q[1];
rz(-pi) q[2];
x q[2];
rz(0.33552334) q[3];
sx q[3];
rz(-2.8918859) q[3];
sx q[3];
rz(1.5229131) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.49332508) q[2];
sx q[2];
rz(-2.2154634) q[2];
sx q[2];
rz(0.43169272) q[2];
rz(-1.4124983) q[3];
sx q[3];
rz(-2.4192211) q[3];
sx q[3];
rz(-0.13599642) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
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
rz(2.7504808) q[0];
sx q[0];
rz(-1.9727805) q[0];
sx q[0];
rz(-3.0294763) q[0];
rz(0.21513367) q[1];
sx q[1];
rz(-1.5810177) q[1];
sx q[1];
rz(-1.1134061) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4556901) q[0];
sx q[0];
rz(-2.6705574) q[0];
sx q[0];
rz(2.5409565) q[0];
x q[1];
rz(-0.31539519) q[2];
sx q[2];
rz(-1.3239667) q[2];
sx q[2];
rz(-2.18404) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.0943162) q[1];
sx q[1];
rz(-0.89135209) q[1];
sx q[1];
rz(1.1761155) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0180842) q[3];
sx q[3];
rz(-1.7637858) q[3];
sx q[3];
rz(2.6393294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.001174288) q[2];
sx q[2];
rz(-2.4089456) q[2];
sx q[2];
rz(-0.12602885) q[2];
rz(-1.0472939) q[3];
sx q[3];
rz(-1.3207366) q[3];
sx q[3];
rz(0.70820156) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4814608) q[0];
sx q[0];
rz(-2.3829057) q[0];
sx q[0];
rz(1.460176) q[0];
rz(1.8966282) q[1];
sx q[1];
rz(-2.0472725) q[1];
sx q[1];
rz(-1.9326899) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9019933) q[0];
sx q[0];
rz(-1.308846) q[0];
sx q[0];
rz(0.082017935) q[0];
x q[1];
rz(-2.2751341) q[2];
sx q[2];
rz(-0.87477113) q[2];
sx q[2];
rz(-0.28878402) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.066684494) q[1];
sx q[1];
rz(-1.9004596) q[1];
sx q[1];
rz(-0.7633022) q[1];
rz(-2.2015757) q[3];
sx q[3];
rz(-1.7490083) q[3];
sx q[3];
rz(2.7959787) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.37844354) q[2];
sx q[2];
rz(-1.903879) q[2];
sx q[2];
rz(1.8219927) q[2];
rz(0.59213263) q[3];
sx q[3];
rz(-1.416128) q[3];
sx q[3];
rz(-3.106451) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2329344) q[0];
sx q[0];
rz(-2.6180551) q[0];
sx q[0];
rz(-1.3611025) q[0];
rz(1.2110442) q[1];
sx q[1];
rz(-2.2369604) q[1];
sx q[1];
rz(-2.7499054) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64205326) q[0];
sx q[0];
rz(-1.5649438) q[0];
sx q[0];
rz(2.7070579) q[0];
rz(-pi) q[1];
rz(-0.45231818) q[2];
sx q[2];
rz(-1.7512133) q[2];
sx q[2];
rz(2.7197321) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.78545953) q[1];
sx q[1];
rz(-1.2411989) q[1];
sx q[1];
rz(2.3258924) q[1];
rz(-pi) q[2];
rz(1.4976981) q[3];
sx q[3];
rz(-2.0953296) q[3];
sx q[3];
rz(0.91367164) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.6212375) q[2];
sx q[2];
rz(-1.7988127) q[2];
sx q[2];
rz(1.3367782) q[2];
rz(2.3796066) q[3];
sx q[3];
rz(-0.31969324) q[3];
sx q[3];
rz(2.3412162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8005463) q[0];
sx q[0];
rz(-0.30277345) q[0];
sx q[0];
rz(2.5706932) q[0];
rz(-1.7123429) q[1];
sx q[1];
rz(-2.0776904) q[1];
sx q[1];
rz(-2.9796519) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.93962651) q[0];
sx q[0];
rz(-1.8585748) q[0];
sx q[0];
rz(-0.94596011) q[0];
rz(-1.5230721) q[2];
sx q[2];
rz(-0.18880162) q[2];
sx q[2];
rz(-1.7392841) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.75441277) q[1];
sx q[1];
rz(-1.3993235) q[1];
sx q[1];
rz(0.13053723) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.52311388) q[3];
sx q[3];
rz(-1.5418846) q[3];
sx q[3];
rz(1.4873963) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.95742115) q[2];
sx q[2];
rz(-2.0364169) q[2];
sx q[2];
rz(-1.7133678) q[2];
rz(1.1994294) q[3];
sx q[3];
rz(-1.0933484) q[3];
sx q[3];
rz(1.7709581) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4409055) q[0];
sx q[0];
rz(-0.4754684) q[0];
sx q[0];
rz(-1.1204002) q[0];
rz(1.3700925) q[1];
sx q[1];
rz(-2.1961828) q[1];
sx q[1];
rz(-0.97074769) q[1];
rz(-2.2489108) q[2];
sx q[2];
rz(-2.953139) q[2];
sx q[2];
rz(2.1072731) q[2];
rz(0.45392848) q[3];
sx q[3];
rz(-0.47149999) q[3];
sx q[3];
rz(-0.75054689) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
