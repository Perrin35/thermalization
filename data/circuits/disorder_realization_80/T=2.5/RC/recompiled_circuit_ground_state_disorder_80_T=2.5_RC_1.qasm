OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.03117938) q[0];
sx q[0];
rz(-0.92209446) q[0];
sx q[0];
rz(2.3715012) q[0];
rz(0.71647477) q[1];
sx q[1];
rz(4.0840277) q[1];
sx q[1];
rz(10.066636) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72015136) q[0];
sx q[0];
rz(-0.23166616) q[0];
sx q[0];
rz(-0.51584824) q[0];
rz(-1.5719169) q[2];
sx q[2];
rz(-1.5693451) q[2];
sx q[2];
rz(0.077479428) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.227046) q[1];
sx q[1];
rz(-1.7265336) q[1];
sx q[1];
rz(-2.7776633) q[1];
rz(-pi) q[2];
rz(-2.1987207) q[3];
sx q[3];
rz(-0.84318765) q[3];
sx q[3];
rz(2.477463) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.1610819) q[2];
sx q[2];
rz(-1.0970205) q[2];
sx q[2];
rz(-0.80992997) q[2];
rz(-0.03446456) q[3];
sx q[3];
rz(-2.4813215) q[3];
sx q[3];
rz(3.0380429) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63475364) q[0];
sx q[0];
rz(-2.6320808) q[0];
sx q[0];
rz(-0.65938812) q[0];
rz(-1.7022645) q[1];
sx q[1];
rz(-1.6292452) q[1];
sx q[1];
rz(-2.6606681) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5372807) q[0];
sx q[0];
rz(-1.3635646) q[0];
sx q[0];
rz(2.0181433) q[0];
x q[1];
rz(2.4056817) q[2];
sx q[2];
rz(-2.1440268) q[2];
sx q[2];
rz(2.5664751) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.7742851) q[1];
sx q[1];
rz(-1.3626422) q[1];
sx q[1];
rz(-0.03117301) q[1];
x q[2];
rz(-0.65851237) q[3];
sx q[3];
rz(-0.69039804) q[3];
sx q[3];
rz(1.8883324) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.3769569) q[2];
sx q[2];
rz(-2.3048293) q[2];
sx q[2];
rz(0.62977201) q[2];
rz(-1.1791641) q[3];
sx q[3];
rz(-2.4383014) q[3];
sx q[3];
rz(1.111697) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69029194) q[0];
sx q[0];
rz(-0.010882219) q[0];
sx q[0];
rz(2.1731398) q[0];
rz(-3.0015266) q[1];
sx q[1];
rz(-1.7882971) q[1];
sx q[1];
rz(-0.5828988) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66630689) q[0];
sx q[0];
rz(-2.9687442) q[0];
sx q[0];
rz(-2.8768507) q[0];
rz(-pi) q[1];
x q[1];
rz(0.83816041) q[2];
sx q[2];
rz(-1.3930905) q[2];
sx q[2];
rz(-3.1237683) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.0576833) q[1];
sx q[1];
rz(-2.1321802) q[1];
sx q[1];
rz(-2.3777005) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8898592) q[3];
sx q[3];
rz(-1.8382065) q[3];
sx q[3];
rz(1.2280994) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.6110903) q[2];
sx q[2];
rz(-0.077433057) q[2];
sx q[2];
rz(0.21406847) q[2];
rz(2.8392082) q[3];
sx q[3];
rz(-0.74367911) q[3];
sx q[3];
rz(-2.7157057) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18836235) q[0];
sx q[0];
rz(-1.0091877) q[0];
sx q[0];
rz(1.0444214) q[0];
rz(-0.87772477) q[1];
sx q[1];
rz(-1.5526155) q[1];
sx q[1];
rz(-2.9728319) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89655748) q[0];
sx q[0];
rz(-1.3447133) q[0];
sx q[0];
rz(-2.079179) q[0];
rz(-pi) q[1];
rz(2.9364768) q[2];
sx q[2];
rz(-1.8913219) q[2];
sx q[2];
rz(-1.4687302) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.61353325) q[1];
sx q[1];
rz(-0.643104) q[1];
sx q[1];
rz(-0.36524857) q[1];
rz(2.2764858) q[3];
sx q[3];
rz(-1.5388515) q[3];
sx q[3];
rz(2.3695996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.3622482) q[2];
sx q[2];
rz(-1.2835953) q[2];
sx q[2];
rz(2.0812422) q[2];
rz(-0.29252163) q[3];
sx q[3];
rz(-0.7258324) q[3];
sx q[3];
rz(-3.0574851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5523858) q[0];
sx q[0];
rz(-2.4942106) q[0];
sx q[0];
rz(-2.2733083) q[0];
rz(0.6262511) q[1];
sx q[1];
rz(-1.3887082) q[1];
sx q[1];
rz(1.7832322) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3247899) q[0];
sx q[0];
rz(-2.8395617) q[0];
sx q[0];
rz(-0.30257757) q[0];
rz(-1.7873853) q[2];
sx q[2];
rz(-1.944245) q[2];
sx q[2];
rz(0.83080705) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-3.061314) q[1];
sx q[1];
rz(-1.3484203) q[1];
sx q[1];
rz(-0.23134065) q[1];
rz(2.5739838) q[3];
sx q[3];
rz(-0.72101147) q[3];
sx q[3];
rz(-2.9182059) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.2061578) q[2];
sx q[2];
rz(-2.3249966) q[2];
sx q[2];
rz(0.52331501) q[2];
rz(0.37007904) q[3];
sx q[3];
rz(-2.3612634) q[3];
sx q[3];
rz(1.3453329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10941457) q[0];
sx q[0];
rz(-2.9258756) q[0];
sx q[0];
rz(-0.35414645) q[0];
rz(-0.94193637) q[1];
sx q[1];
rz(-1.6862005) q[1];
sx q[1];
rz(1.8249493) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38625941) q[0];
sx q[0];
rz(-2.0830031) q[0];
sx q[0];
rz(-1.451158) q[0];
rz(-0.022158547) q[2];
sx q[2];
rz(-0.54931927) q[2];
sx q[2];
rz(-1.4469128) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.9256014) q[1];
sx q[1];
rz(-0.18454889) q[1];
sx q[1];
rz(2.8021003) q[1];
x q[2];
rz(0.64423465) q[3];
sx q[3];
rz(-1.3431864) q[3];
sx q[3];
rz(-1.5952974) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.3256623) q[2];
sx q[2];
rz(-0.38620913) q[2];
sx q[2];
rz(-0.33829921) q[2];
rz(-2.6650688) q[3];
sx q[3];
rz(-0.75348133) q[3];
sx q[3];
rz(2.8009955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7020096) q[0];
sx q[0];
rz(-1.5583353) q[0];
sx q[0];
rz(-0.53771341) q[0];
rz(-1.733755) q[1];
sx q[1];
rz(-2.6815963) q[1];
sx q[1];
rz(-2.5163311) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0940373) q[0];
sx q[0];
rz(-0.68935822) q[0];
sx q[0];
rz(-2.5563142) q[0];
rz(1.6444667) q[2];
sx q[2];
rz(-2.1084614) q[2];
sx q[2];
rz(1.0755838) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.3846674) q[1];
sx q[1];
rz(-2.1672492) q[1];
sx q[1];
rz(1.3104964) q[1];
x q[2];
rz(2.7688249) q[3];
sx q[3];
rz(-2.3298752) q[3];
sx q[3];
rz(1.2889287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.9739146) q[2];
sx q[2];
rz(-0.46678552) q[2];
sx q[2];
rz(-1.3803049) q[2];
rz(0.4365094) q[3];
sx q[3];
rz(-2.1197539) q[3];
sx q[3];
rz(0.71353394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1768782) q[0];
sx q[0];
rz(-0.62129337) q[0];
sx q[0];
rz(-0.51625133) q[0];
rz(2.3618354) q[1];
sx q[1];
rz(-0.98194352) q[1];
sx q[1];
rz(1.0468743) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4488432) q[0];
sx q[0];
rz(-1.7147281) q[0];
sx q[0];
rz(2.8400665) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4372836) q[2];
sx q[2];
rz(-0.71137911) q[2];
sx q[2];
rz(-2.7976183) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.88491625) q[1];
sx q[1];
rz(-3.0312928) q[1];
sx q[1];
rz(-2.1259456) q[1];
rz(-2.6627296) q[3];
sx q[3];
rz(-1.2148464) q[3];
sx q[3];
rz(1.009481) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.20389916) q[2];
sx q[2];
rz(-0.1940618) q[2];
sx q[2];
rz(-0.95721179) q[2];
rz(0.29414487) q[3];
sx q[3];
rz(-1.1295986) q[3];
sx q[3];
rz(0.2778151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
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
rz(0.99943632) q[0];
sx q[0];
rz(-2.9111828) q[0];
sx q[0];
rz(-2.9545422) q[0];
rz(-2.710178) q[1];
sx q[1];
rz(-2.7038733) q[1];
sx q[1];
rz(-1.6492708) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6607912) q[0];
sx q[0];
rz(-1.6422762) q[0];
sx q[0];
rz(0.070822318) q[0];
rz(1.0643105) q[2];
sx q[2];
rz(-0.37831719) q[2];
sx q[2];
rz(-1.0418721) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.3826661) q[1];
sx q[1];
rz(-1.1622283) q[1];
sx q[1];
rz(2.9620785) q[1];
rz(-pi) q[2];
x q[2];
rz(0.093676693) q[3];
sx q[3];
rz(-1.7266577) q[3];
sx q[3];
rz(2.3222498) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6734068) q[2];
sx q[2];
rz(-0.91172051) q[2];
sx q[2];
rz(2.1569596) q[2];
rz(-0.51472384) q[3];
sx q[3];
rz(-2.616021) q[3];
sx q[3];
rz(-1.908186) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.89727) q[0];
sx q[0];
rz(-1.4800625) q[0];
sx q[0];
rz(0.75862128) q[0];
rz(1.9482535) q[1];
sx q[1];
rz(-1.1921644) q[1];
sx q[1];
rz(1.4512482) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39501909) q[0];
sx q[0];
rz(-1.7749471) q[0];
sx q[0];
rz(1.3373119) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.042165857) q[2];
sx q[2];
rz(-1.773196) q[2];
sx q[2];
rz(-0.037994904) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9368694) q[1];
sx q[1];
rz(-1.544049) q[1];
sx q[1];
rz(-0.95519798) q[1];
rz(-2.8192807) q[3];
sx q[3];
rz(-1.1813643) q[3];
sx q[3];
rz(-0.001231391) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.5997368) q[2];
sx q[2];
rz(-2.2208322) q[2];
sx q[2];
rz(2.3403781) q[2];
rz(-2.2090705) q[3];
sx q[3];
rz(-1.2152117) q[3];
sx q[3];
rz(2.7264989) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[3];
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
rz(1.5634609) q[0];
sx q[0];
rz(-1.3787855) q[0];
sx q[0];
rz(-0.61510573) q[0];
rz(0.32147944) q[1];
sx q[1];
rz(-2.165806) q[1];
sx q[1];
rz(1.4478366) q[1];
rz(1.4463439) q[2];
sx q[2];
rz(-1.4067408) q[2];
sx q[2];
rz(-0.066700145) q[2];
rz(0.23033167) q[3];
sx q[3];
rz(-1.2970222) q[3];
sx q[3];
rz(2.6877689) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
