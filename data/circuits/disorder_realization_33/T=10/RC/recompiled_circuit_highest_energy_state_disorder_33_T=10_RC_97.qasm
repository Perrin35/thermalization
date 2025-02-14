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
rz(-3.1193982) q[0];
sx q[0];
rz(-0.57214195) q[0];
sx q[0];
rz(-0.194508) q[0];
rz(-1.0021078) q[1];
sx q[1];
rz(-0.78295308) q[1];
sx q[1];
rz(-0.015838239) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7730071) q[0];
sx q[0];
rz(-1.9814648) q[0];
sx q[0];
rz(2.3334614) q[0];
x q[1];
rz(1.1646526) q[2];
sx q[2];
rz(-2.32562) q[2];
sx q[2];
rz(1.3938122) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.046151) q[1];
sx q[1];
rz(-2.4956551) q[1];
sx q[1];
rz(-1.7424042) q[1];
x q[2];
rz(1.2239816) q[3];
sx q[3];
rz(-1.6458316) q[3];
sx q[3];
rz(-1.1574367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.013449239) q[2];
sx q[2];
rz(-1.6703419) q[2];
sx q[2];
rz(0.61689287) q[2];
rz(-1.8938176) q[3];
sx q[3];
rz(-2.3168677) q[3];
sx q[3];
rz(2.8491546) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1013252) q[0];
sx q[0];
rz(-1.9693002) q[0];
sx q[0];
rz(-0.030666703) q[0];
rz(-1.4847633) q[1];
sx q[1];
rz(-2.8430884) q[1];
sx q[1];
rz(-0.31918496) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6761665) q[0];
sx q[0];
rz(-0.54252189) q[0];
sx q[0];
rz(-1.2529621) q[0];
rz(1.9640552) q[2];
sx q[2];
rz(-1.3280091) q[2];
sx q[2];
rz(0.78906203) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.044925) q[1];
sx q[1];
rz(-0.54122347) q[1];
sx q[1];
rz(0.9344395) q[1];
x q[2];
rz(-2.3461211) q[3];
sx q[3];
rz(-0.69517259) q[3];
sx q[3];
rz(0.89423394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.9926051) q[2];
sx q[2];
rz(-0.4036029) q[2];
sx q[2];
rz(2.4888424) q[2];
rz(2.2720689) q[3];
sx q[3];
rz(-1.0966938) q[3];
sx q[3];
rz(-2.7147527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2000548) q[0];
sx q[0];
rz(-2.8917199) q[0];
sx q[0];
rz(0.13724929) q[0];
rz(-0.50357729) q[1];
sx q[1];
rz(-0.63882393) q[1];
sx q[1];
rz(0.30127475) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6414911) q[0];
sx q[0];
rz(-1.8805686) q[0];
sx q[0];
rz(-0.70112438) q[0];
rz(1.7687479) q[2];
sx q[2];
rz(-2.8471409) q[2];
sx q[2];
rz(2.1535923) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.2259146) q[1];
sx q[1];
rz(-1.3667733) q[1];
sx q[1];
rz(-1.4218487) q[1];
rz(-0.6676612) q[3];
sx q[3];
rz(-1.6958837) q[3];
sx q[3];
rz(1.1043105) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0425903) q[2];
sx q[2];
rz(-0.66163915) q[2];
sx q[2];
rz(-2.674687) q[2];
rz(0.0070988797) q[3];
sx q[3];
rz(-0.10741281) q[3];
sx q[3];
rz(-0.98889178) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0840833) q[0];
sx q[0];
rz(-3.1355661) q[0];
sx q[0];
rz(-0.65473336) q[0];
rz(-1.5714802) q[1];
sx q[1];
rz(-1.7127769) q[1];
sx q[1];
rz(-0.44081259) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8705298) q[0];
sx q[0];
rz(-1.5952627) q[0];
sx q[0];
rz(-1.5951675) q[0];
rz(0.019182656) q[2];
sx q[2];
rz(-0.47878107) q[2];
sx q[2];
rz(-2.5342902) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.1487351) q[1];
sx q[1];
rz(-1.8655708) q[1];
sx q[1];
rz(1.8026505) q[1];
rz(-pi) q[2];
rz(1.8809404) q[3];
sx q[3];
rz(-2.0278483) q[3];
sx q[3];
rz(3.0503805) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.4948027) q[2];
sx q[2];
rz(-1.9014634) q[2];
sx q[2];
rz(-0.49171641) q[2];
rz(-1.3382737) q[3];
sx q[3];
rz(-1.0889784) q[3];
sx q[3];
rz(1.6751539) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3717475) q[0];
sx q[0];
rz(-1.0862229) q[0];
sx q[0];
rz(1.0787429) q[0];
rz(3.0951913) q[1];
sx q[1];
rz(-1.6175783) q[1];
sx q[1];
rz(-2.0414415) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85202269) q[0];
sx q[0];
rz(-1.6873813) q[0];
sx q[0];
rz(1.5615211) q[0];
x q[1];
rz(-0.64952594) q[2];
sx q[2];
rz(-2.1500271) q[2];
sx q[2];
rz(-1.1729826) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8535159) q[1];
sx q[1];
rz(-1.8156681) q[1];
sx q[1];
rz(0.71611407) q[1];
rz(-pi) q[2];
rz(-1.4579822) q[3];
sx q[3];
rz(-1.3028909) q[3];
sx q[3];
rz(-2.9873141) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.2940353) q[2];
sx q[2];
rz(-0.39176771) q[2];
sx q[2];
rz(0.34226391) q[2];
rz(1.3058454) q[3];
sx q[3];
rz(-1.9394453) q[3];
sx q[3];
rz(2.0391298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-0.24868988) q[0];
sx q[0];
rz(-1.1868287) q[0];
sx q[0];
rz(2.7500395) q[0];
rz(-0.64741099) q[1];
sx q[1];
rz(-2.7029111) q[1];
sx q[1];
rz(-0.90519261) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0150512) q[0];
sx q[0];
rz(-0.76989664) q[0];
sx q[0];
rz(-0.73791821) q[0];
rz(-pi) q[1];
rz(-0.91475418) q[2];
sx q[2];
rz(-1.1335271) q[2];
sx q[2];
rz(0.8503051) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1282922) q[1];
sx q[1];
rz(-0.97620539) q[1];
sx q[1];
rz(1.250729) q[1];
x q[2];
rz(-2.9764683) q[3];
sx q[3];
rz(-2.8493463) q[3];
sx q[3];
rz(-1.0632893) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.68760005) q[2];
sx q[2];
rz(-0.31012055) q[2];
sx q[2];
rz(1.0529168) q[2];
rz(2.8478801) q[3];
sx q[3];
rz(-0.83455825) q[3];
sx q[3];
rz(2.6101051) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76327908) q[0];
sx q[0];
rz(-1.8623619) q[0];
sx q[0];
rz(-3.1148425) q[0];
rz(-1.1862952) q[1];
sx q[1];
rz(-0.18272884) q[1];
sx q[1];
rz(0.15348405) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5993921) q[0];
sx q[0];
rz(-1.8398823) q[0];
sx q[0];
rz(-2.0890636) q[0];
rz(-pi) q[1];
rz(-0.085603733) q[2];
sx q[2];
rz(-0.55194469) q[2];
sx q[2];
rz(-2.7526223) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3378422) q[1];
sx q[1];
rz(-2.1598247) q[1];
sx q[1];
rz(-2.4250125) q[1];
rz(-pi) q[2];
x q[2];
rz(0.028702486) q[3];
sx q[3];
rz(-1.5099635) q[3];
sx q[3];
rz(-2.3059152) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.2457876) q[2];
sx q[2];
rz(-0.20496002) q[2];
sx q[2];
rz(1.2650355) q[2];
rz(-2.9686019) q[3];
sx q[3];
rz(-1.750662) q[3];
sx q[3];
rz(-1.0783819) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
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
rz(-3.0377401) q[0];
sx q[0];
rz(-0.2094035) q[0];
sx q[0];
rz(-3.0331392) q[0];
rz(1.3718038) q[1];
sx q[1];
rz(-1.7864952) q[1];
sx q[1];
rz(-3.1350737) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9193503) q[0];
sx q[0];
rz(-2.0323638) q[0];
sx q[0];
rz(0.95141769) q[0];
rz(-1.4638607) q[2];
sx q[2];
rz(-1.7902014) q[2];
sx q[2];
rz(2.1768513) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.5287787) q[1];
sx q[1];
rz(-1.0135796) q[1];
sx q[1];
rz(0.91220906) q[1];
rz(1.9912854) q[3];
sx q[3];
rz(-0.83339993) q[3];
sx q[3];
rz(-2.5954136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6846377) q[2];
sx q[2];
rz(-2.4943887) q[2];
sx q[2];
rz(1.7015438) q[2];
rz(2.7889934) q[3];
sx q[3];
rz(-2.2779901) q[3];
sx q[3];
rz(0.45261228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0532992) q[0];
sx q[0];
rz(-2.0501917) q[0];
sx q[0];
rz(-1.1540867) q[0];
rz(-1.6905009) q[1];
sx q[1];
rz(-2.4595551) q[1];
sx q[1];
rz(2.3114204) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51541942) q[0];
sx q[0];
rz(-0.011224482) q[0];
sx q[0];
rz(-2.8521538) q[0];
x q[1];
rz(-2.9014456) q[2];
sx q[2];
rz(-2.6278751) q[2];
sx q[2];
rz(1.7983537) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.3026003) q[1];
sx q[1];
rz(-2.5328089) q[1];
sx q[1];
rz(-2.3997612) q[1];
rz(-pi) q[2];
rz(0.19751246) q[3];
sx q[3];
rz(-2.2220051) q[3];
sx q[3];
rz(0.75266389) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.1900968) q[2];
sx q[2];
rz(-2.2647965) q[2];
sx q[2];
rz(-2.4693176) q[2];
rz(-2.7049474) q[3];
sx q[3];
rz(-1.838622) q[3];
sx q[3];
rz(0.28948998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.745568) q[0];
sx q[0];
rz(-0.73306274) q[0];
sx q[0];
rz(-0.42098862) q[0];
rz(1.9898532) q[1];
sx q[1];
rz(-0.20944171) q[1];
sx q[1];
rz(-2.4236325) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5508681) q[0];
sx q[0];
rz(-1.8483254) q[0];
sx q[0];
rz(3.1184349) q[0];
x q[1];
rz(0.69952632) q[2];
sx q[2];
rz(-1.8803863) q[2];
sx q[2];
rz(1.6242611) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.5567498) q[1];
sx q[1];
rz(-1.6989917) q[1];
sx q[1];
rz(1.8845647) q[1];
rz(2.8676492) q[3];
sx q[3];
rz(-2.0060325) q[3];
sx q[3];
rz(2.1147605) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(3.1154321) q[2];
sx q[2];
rz(-0.60439503) q[2];
sx q[2];
rz(-1.3755414) q[2];
rz(-0.17624217) q[3];
sx q[3];
rz(-2.3459489) q[3];
sx q[3];
rz(3.0680883) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7701223) q[0];
sx q[0];
rz(-1.4513411) q[0];
sx q[0];
rz(2.0392192) q[0];
rz(2.2813028) q[1];
sx q[1];
rz(-1.6382891) q[1];
sx q[1];
rz(2.0232497) q[1];
rz(-0.76535759) q[2];
sx q[2];
rz(-1.6156964) q[2];
sx q[2];
rz(2.2888301) q[2];
rz(2.7793814) q[3];
sx q[3];
rz(-1.9345761) q[3];
sx q[3];
rz(-3.0349845) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
