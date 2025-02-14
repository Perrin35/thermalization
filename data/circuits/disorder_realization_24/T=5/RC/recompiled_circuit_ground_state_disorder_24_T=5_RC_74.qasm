OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.71819031) q[0];
sx q[0];
rz(-3.1388404) q[0];
sx q[0];
rz(1.0116853) q[0];
rz(-2.6818795) q[1];
sx q[1];
rz(-1.3016394) q[1];
sx q[1];
rz(0.17170061) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26352208) q[0];
sx q[0];
rz(-1.4891013) q[0];
sx q[0];
rz(-0.6058713) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6284211) q[2];
sx q[2];
rz(-2.9193239) q[2];
sx q[2];
rz(1.5813511) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.2312647) q[1];
sx q[1];
rz(-0.99580169) q[1];
sx q[1];
rz(2.0851589) q[1];
rz(1.4775253) q[3];
sx q[3];
rz(-1.516946) q[3];
sx q[3];
rz(-0.71291956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1823938) q[2];
sx q[2];
rz(-0.48754075) q[2];
sx q[2];
rz(-1.7215151) q[2];
rz(-1.9567664) q[3];
sx q[3];
rz(-1.5337557) q[3];
sx q[3];
rz(-2.0737295) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
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
rz(-0.083953388) q[0];
sx q[0];
rz(-1.5204484) q[0];
sx q[0];
rz(0.49348304) q[0];
rz(2.3656942) q[1];
sx q[1];
rz(-0.50965613) q[1];
sx q[1];
rz(2.6367771) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0836661) q[0];
sx q[0];
rz(-0.89910451) q[0];
sx q[0];
rz(-0.85606411) q[0];
rz(-pi) q[1];
rz(2.6569064) q[2];
sx q[2];
rz(-0.58161608) q[2];
sx q[2];
rz(1.5145375) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.3146554) q[1];
sx q[1];
rz(-0.75410226) q[1];
sx q[1];
rz(-2.5376863) q[1];
x q[2];
rz(2.9166031) q[3];
sx q[3];
rz(-0.27249042) q[3];
sx q[3];
rz(-0.30411965) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.8043171) q[2];
sx q[2];
rz(-1.8307468) q[2];
sx q[2];
rz(-0.47067019) q[2];
rz(-0.63052952) q[3];
sx q[3];
rz(-2.1157406) q[3];
sx q[3];
rz(1.7414909) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.077496342) q[0];
sx q[0];
rz(-0.12501669) q[0];
sx q[0];
rz(-0.65573829) q[0];
rz(-2.8302622) q[1];
sx q[1];
rz(-0.89404023) q[1];
sx q[1];
rz(0.071455926) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27256672) q[0];
sx q[0];
rz(-1.9211968) q[0];
sx q[0];
rz(-3.0883771) q[0];
x q[1];
rz(-1.2024443) q[2];
sx q[2];
rz(-1.2542274) q[2];
sx q[2];
rz(-0.20395111) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.25117043) q[1];
sx q[1];
rz(-2.648472) q[1];
sx q[1];
rz(-1.432073) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7830332) q[3];
sx q[3];
rz(-2.0925412) q[3];
sx q[3];
rz(2.9468342) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.91325703) q[2];
sx q[2];
rz(-1.5847289) q[2];
sx q[2];
rz(-0.060001686) q[2];
rz(-1.9289198) q[3];
sx q[3];
rz(-0.69827497) q[3];
sx q[3];
rz(-0.76446271) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-0.54670984) q[0];
sx q[0];
rz(-1.3285652) q[0];
sx q[0];
rz(0.57719624) q[0];
rz(-2.3120841) q[1];
sx q[1];
rz(-1.0521768) q[1];
sx q[1];
rz(-0.83782354) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4234377) q[0];
sx q[0];
rz(-2.1912247) q[0];
sx q[0];
rz(1.4761094) q[0];
x q[1];
rz(-2.1609862) q[2];
sx q[2];
rz(-1.6211121) q[2];
sx q[2];
rz(-0.8129763) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.438617) q[1];
sx q[1];
rz(-1.9432406) q[1];
sx q[1];
rz(2.245642) q[1];
rz(-pi) q[2];
rz(-0.44562037) q[3];
sx q[3];
rz(-0.94404781) q[3];
sx q[3];
rz(2.8901951) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.025658) q[2];
sx q[2];
rz(-1.5237153) q[2];
sx q[2];
rz(-1.9177829) q[2];
rz(0.4661679) q[3];
sx q[3];
rz(-0.40188447) q[3];
sx q[3];
rz(0.60552067) q[3];
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
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25662988) q[0];
sx q[0];
rz(-1.7054568) q[0];
sx q[0];
rz(-0.10109854) q[0];
rz(-0.413232) q[1];
sx q[1];
rz(-2.7006472) q[1];
sx q[1];
rz(-1.0879999) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1032858) q[0];
sx q[0];
rz(-1.8977471) q[0];
sx q[0];
rz(1.0590382) q[0];
rz(-pi) q[1];
rz(-0.83482439) q[2];
sx q[2];
rz(-1.2918279) q[2];
sx q[2];
rz(0.5912515) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.99876304) q[1];
sx q[1];
rz(-1.8149477) q[1];
sx q[1];
rz(-1.2296089) q[1];
rz(-pi) q[2];
rz(-1.5002046) q[3];
sx q[3];
rz(-0.90259457) q[3];
sx q[3];
rz(-1.6833504) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.2130412) q[2];
sx q[2];
rz(-0.88625675) q[2];
sx q[2];
rz(1.586033) q[2];
rz(-2.0689615) q[3];
sx q[3];
rz(-1.5242679) q[3];
sx q[3];
rz(-0.13256375) q[3];
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
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9685386) q[0];
sx q[0];
rz(-2.2569077) q[0];
sx q[0];
rz(-1.7065077) q[0];
rz(-0.75812078) q[1];
sx q[1];
rz(-0.70223141) q[1];
sx q[1];
rz(-0.10163669) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.223343) q[0];
sx q[0];
rz(-2.0416284) q[0];
sx q[0];
rz(0.74145326) q[0];
rz(2.7804537) q[2];
sx q[2];
rz(-1.9996579) q[2];
sx q[2];
rz(2.7782235) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.95127997) q[1];
sx q[1];
rz(-1.4644831) q[1];
sx q[1];
rz(-2.8578958) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6061344) q[3];
sx q[3];
rz(-2.1369918) q[3];
sx q[3];
rz(-2.9009852) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.8018084) q[2];
sx q[2];
rz(-1.1016176) q[2];
sx q[2];
rz(1.3327117) q[2];
rz(-2.0594275) q[3];
sx q[3];
rz(-0.75282955) q[3];
sx q[3];
rz(1.4235206) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4448755) q[0];
sx q[0];
rz(-1.8901905) q[0];
sx q[0];
rz(0.74309293) q[0];
rz(-2.3785059) q[1];
sx q[1];
rz(-0.63947314) q[1];
sx q[1];
rz(0.9185763) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0325286) q[0];
sx q[0];
rz(-1.5957513) q[0];
sx q[0];
rz(0.0080913261) q[0];
x q[1];
rz(0.95942504) q[2];
sx q[2];
rz(-0.55771962) q[2];
sx q[2];
rz(-3.0447247) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.1700701) q[1];
sx q[1];
rz(-0.60107175) q[1];
sx q[1];
rz(-2.5843589) q[1];
x q[2];
rz(3.0094163) q[3];
sx q[3];
rz(-2.4130947) q[3];
sx q[3];
rz(0.26154172) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.7243728) q[2];
sx q[2];
rz(-1.1823187) q[2];
sx q[2];
rz(2.7057538) q[2];
rz(3.0271652) q[3];
sx q[3];
rz(-2.8947713) q[3];
sx q[3];
rz(0.59823263) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8082387) q[0];
sx q[0];
rz(-0.55753189) q[0];
sx q[0];
rz(0.58407855) q[0];
rz(0.43560478) q[1];
sx q[1];
rz(-1.4608811) q[1];
sx q[1];
rz(1.6216507) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7652033) q[0];
sx q[0];
rz(-1.1113154) q[0];
sx q[0];
rz(-2.8110519) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9797835) q[2];
sx q[2];
rz(-1.7715381) q[2];
sx q[2];
rz(0.44524064) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.4055109) q[1];
sx q[1];
rz(-0.20624948) q[1];
sx q[1];
rz(-3.095605) q[1];
rz(2.5831843) q[3];
sx q[3];
rz(-0.30277751) q[3];
sx q[3];
rz(-1.4608135) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(3.0674151) q[2];
sx q[2];
rz(-1.7895074) q[2];
sx q[2];
rz(1.1161067) q[2];
rz(1.3793147) q[3];
sx q[3];
rz(-1.1032871) q[3];
sx q[3];
rz(-0.11837676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8922358) q[0];
sx q[0];
rz(-3.0638969) q[0];
sx q[0];
rz(-0.78999162) q[0];
rz(-0.063591592) q[1];
sx q[1];
rz(-0.60706943) q[1];
sx q[1];
rz(2.7008609) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1871619) q[0];
sx q[0];
rz(-1.2557898) q[0];
sx q[0];
rz(-0.31401547) q[0];
rz(-pi) q[1];
x q[1];
rz(0.49479167) q[2];
sx q[2];
rz(-1.2107009) q[2];
sx q[2];
rz(0.39967135) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.4197582) q[1];
sx q[1];
rz(-1.7350619) q[1];
sx q[1];
rz(2.052015) q[1];
rz(-pi) q[2];
x q[2];
rz(0.041429511) q[3];
sx q[3];
rz(-1.7600312) q[3];
sx q[3];
rz(1.7570868) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.26238394) q[2];
sx q[2];
rz(-2.3449506) q[2];
sx q[2];
rz(0.66514307) q[2];
rz(2.3696259) q[3];
sx q[3];
rz(-0.77163458) q[3];
sx q[3];
rz(1.6316679) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36122286) q[0];
sx q[0];
rz(-0.64798111) q[0];
sx q[0];
rz(1.004647) q[0];
rz(1.4077582) q[1];
sx q[1];
rz(-1.6561457) q[1];
sx q[1];
rz(-1.2900603) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4125975) q[0];
sx q[0];
rz(-2.3069803) q[0];
sx q[0];
rz(2.8060629) q[0];
rz(2.9384841) q[2];
sx q[2];
rz(-0.98047335) q[2];
sx q[2];
rz(-0.029435722) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.87566808) q[1];
sx q[1];
rz(-2.7920715) q[1];
sx q[1];
rz(-2.4825579) q[1];
x q[2];
rz(-2.9674888) q[3];
sx q[3];
rz(-2.1645567) q[3];
sx q[3];
rz(-1.6792149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.1349692) q[2];
sx q[2];
rz(-1.1639736) q[2];
sx q[2];
rz(0.25035614) q[2];
rz(0.97744673) q[3];
sx q[3];
rz(-1.9340632) q[3];
sx q[3];
rz(-1.6710963) q[3];
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
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95494315) q[0];
sx q[0];
rz(-1.432812) q[0];
sx q[0];
rz(-1.1695255) q[0];
rz(2.3969338) q[1];
sx q[1];
rz(-2.3458993) q[1];
sx q[1];
rz(0.69534272) q[1];
rz(1.5834687) q[2];
sx q[2];
rz(-1.7141533) q[2];
sx q[2];
rz(-1.2319596) q[2];
rz(-0.43375265) q[3];
sx q[3];
rz(-2.0278145) q[3];
sx q[3];
rz(-0.47587004) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
