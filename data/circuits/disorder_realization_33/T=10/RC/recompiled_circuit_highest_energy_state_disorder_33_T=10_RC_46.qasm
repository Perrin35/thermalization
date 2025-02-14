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
rz(2.1394849) q[1];
sx q[1];
rz(-2.3586396) q[1];
sx q[1];
rz(0.015838239) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.16222157) q[0];
sx q[0];
rz(-0.88484513) q[0];
sx q[0];
rz(2.5995194) q[0];
x q[1];
rz(2.7439666) q[2];
sx q[2];
rz(-0.83774746) q[2];
sx q[2];
rz(0.83329569) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.33783198) q[1];
sx q[1];
rz(-1.4678218) q[1];
sx q[1];
rz(0.93195685) q[1];
x q[2];
rz(0.079766131) q[3];
sx q[3];
rz(-1.9165943) q[3];
sx q[3];
rz(-2.7553204) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.013449239) q[2];
sx q[2];
rz(-1.6703419) q[2];
sx q[2];
rz(0.61689287) q[2];
rz(1.2477751) q[3];
sx q[3];
rz(-0.82472491) q[3];
sx q[3];
rz(0.29243803) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(3.1013252) q[0];
sx q[0];
rz(-1.9693002) q[0];
sx q[0];
rz(-3.110926) q[0];
rz(-1.6568294) q[1];
sx q[1];
rz(-2.8430884) q[1];
sx q[1];
rz(-2.8224077) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3108513) q[0];
sx q[0];
rz(-1.7328528) q[0];
sx q[0];
rz(1.0507163) q[0];
x q[1];
rz(0.99700751) q[2];
sx q[2];
rz(-2.6827723) q[2];
sx q[2];
rz(-0.25653893) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-3.1029237) q[1];
sx q[1];
rz(-1.2596411) q[1];
sx q[1];
rz(1.1204847) q[1];
rz(-pi) q[2];
rz(2.3461211) q[3];
sx q[3];
rz(-0.69517259) q[3];
sx q[3];
rz(-0.89423394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.14898758) q[2];
sx q[2];
rz(-0.4036029) q[2];
sx q[2];
rz(2.4888424) q[2];
rz(-2.2720689) q[3];
sx q[3];
rz(-2.0448989) q[3];
sx q[3];
rz(-2.7147527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.94153786) q[0];
sx q[0];
rz(-0.24987276) q[0];
sx q[0];
rz(0.13724929) q[0];
rz(-0.50357729) q[1];
sx q[1];
rz(-2.5027687) q[1];
sx q[1];
rz(2.8403179) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27583081) q[0];
sx q[0];
rz(-2.3858285) q[0];
sx q[0];
rz(-2.6810032) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7687479) q[2];
sx q[2];
rz(-2.8471409) q[2];
sx q[2];
rz(-2.1535923) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.5884389) q[1];
sx q[1];
rz(-2.8895973) q[1];
sx q[1];
rz(2.5193922) q[1];
x q[2];
rz(-0.6676612) q[3];
sx q[3];
rz(-1.4457089) q[3];
sx q[3];
rz(2.0372822) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.0425903) q[2];
sx q[2];
rz(-2.4799535) q[2];
sx q[2];
rz(0.46690565) q[2];
rz(3.1344938) q[3];
sx q[3];
rz(-3.0341798) q[3];
sx q[3];
rz(-0.98889178) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0840833) q[0];
sx q[0];
rz(-3.1355661) q[0];
sx q[0];
rz(-2.4868593) q[0];
rz(-1.5714802) q[1];
sx q[1];
rz(-1.4288158) q[1];
sx q[1];
rz(0.44081259) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0546608) q[0];
sx q[0];
rz(-3.107061) q[0];
sx q[0];
rz(-2.3582929) q[0];
rz(-pi) q[1];
rz(1.5608403) q[2];
sx q[2];
rz(-1.0921108) q[2];
sx q[2];
rz(-2.5559024) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7880169) q[1];
sx q[1];
rz(-1.792479) q[1];
sx q[1];
rz(-2.8391929) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2606523) q[3];
sx q[3];
rz(-2.0278483) q[3];
sx q[3];
rz(-3.0503805) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.4948027) q[2];
sx q[2];
rz(-1.9014634) q[2];
sx q[2];
rz(2.6498762) q[2];
rz(1.3382737) q[3];
sx q[3];
rz(-1.0889784) q[3];
sx q[3];
rz(-1.6751539) q[3];
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
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76984513) q[0];
sx q[0];
rz(-1.0862229) q[0];
sx q[0];
rz(-1.0787429) q[0];
rz(-0.046401333) q[1];
sx q[1];
rz(-1.6175783) q[1];
sx q[1];
rz(-2.0414415) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3691424) q[0];
sx q[0];
rz(-3.0246409) q[0];
sx q[0];
rz(3.0625615) q[0];
x q[1];
rz(-2.3170839) q[2];
sx q[2];
rz(-2.3002599) q[2];
sx q[2];
rz(2.1192604) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-3.0667923) q[1];
sx q[1];
rz(-0.88035816) q[1];
sx q[1];
rz(1.2509173) q[1];
rz(-pi) q[2];
rz(0.2695378) q[3];
sx q[3];
rz(-1.4620228) q[3];
sx q[3];
rz(-1.6950932) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.84755737) q[2];
sx q[2];
rz(-0.39176771) q[2];
sx q[2];
rz(-2.7993287) q[2];
rz(1.8357473) q[3];
sx q[3];
rz(-1.9394453) q[3];
sx q[3];
rz(-2.0391298) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24868988) q[0];
sx q[0];
rz(-1.9547639) q[0];
sx q[0];
rz(2.7500395) q[0];
rz(2.4941817) q[1];
sx q[1];
rz(-2.7029111) q[1];
sx q[1];
rz(2.2364) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0150512) q[0];
sx q[0];
rz(-2.371696) q[0];
sx q[0];
rz(0.73791821) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6086285) q[2];
sx q[2];
rz(-0.98532444) q[2];
sx q[2];
rz(1.0356071) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.47899095) q[1];
sx q[1];
rz(-0.66598611) q[1];
sx q[1];
rz(-2.7061092) q[1];
rz(-pi) q[2];
rz(-2.9764683) q[3];
sx q[3];
rz(-0.29224631) q[3];
sx q[3];
rz(1.0632893) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.4539926) q[2];
sx q[2];
rz(-0.31012055) q[2];
sx q[2];
rz(-1.0529168) q[2];
rz(-2.8478801) q[3];
sx q[3];
rz(-2.3070344) q[3];
sx q[3];
rz(2.6101051) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76327908) q[0];
sx q[0];
rz(-1.2792307) q[0];
sx q[0];
rz(-0.026750201) q[0];
rz(-1.1862952) q[1];
sx q[1];
rz(-2.9588638) q[1];
sx q[1];
rz(-0.15348405) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1790548) q[0];
sx q[0];
rz(-2.0686596) q[0];
sx q[0];
rz(0.30740096) q[0];
rz(3.0559889) q[2];
sx q[2];
rz(-2.589648) q[2];
sx q[2];
rz(-0.38897038) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.9423034) q[1];
sx q[1];
rz(-0.89312141) q[1];
sx q[1];
rz(-0.79395644) q[1];
rz(-pi) q[2];
rz(-3.1128902) q[3];
sx q[3];
rz(-1.5099635) q[3];
sx q[3];
rz(0.83567747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.89580506) q[2];
sx q[2];
rz(-0.20496002) q[2];
sx q[2];
rz(1.2650355) q[2];
rz(0.17299077) q[3];
sx q[3];
rz(-1.3909307) q[3];
sx q[3];
rz(1.0783819) q[3];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10385253) q[0];
sx q[0];
rz(-0.2094035) q[0];
sx q[0];
rz(0.10845342) q[0];
rz(1.7697889) q[1];
sx q[1];
rz(-1.7864952) q[1];
sx q[1];
rz(-0.0065189204) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9193503) q[0];
sx q[0];
rz(-2.0323638) q[0];
sx q[0];
rz(2.190175) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.22062515) q[2];
sx q[2];
rz(-1.675159) q[2];
sx q[2];
rz(-2.5121784) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7951736) q[1];
sx q[1];
rz(-1.0246207) q[1];
sx q[1];
rz(0.66728981) q[1];
rz(-pi) q[2];
rz(1.9912854) q[3];
sx q[3];
rz(-0.83339993) q[3];
sx q[3];
rz(-2.5954136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.6846377) q[2];
sx q[2];
rz(-0.64720398) q[2];
sx q[2];
rz(1.4400488) q[2];
rz(2.7889934) q[3];
sx q[3];
rz(-0.86360258) q[3];
sx q[3];
rz(2.6889804) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.088293485) q[0];
sx q[0];
rz(-2.0501917) q[0];
sx q[0];
rz(1.1540867) q[0];
rz(-1.6905009) q[1];
sx q[1];
rz(-0.68203753) q[1];
sx q[1];
rz(0.83017224) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.80487554) q[0];
sx q[0];
rz(-1.5600388) q[0];
sx q[0];
rz(-1.5740001) q[0];
rz(-pi) q[1];
rz(-2.9014456) q[2];
sx q[2];
rz(-0.51371759) q[2];
sx q[2];
rz(-1.7983537) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.37646078) q[1];
sx q[1];
rz(-1.174095) q[1];
sx q[1];
rz(-2.6668616) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8228047) q[3];
sx q[3];
rz(-2.4652836) q[3];
sx q[3];
rz(-2.0700434) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.95149583) q[2];
sx q[2];
rz(-2.2647965) q[2];
sx q[2];
rz(0.67227501) q[2];
rz(0.43664524) q[3];
sx q[3];
rz(-1.838622) q[3];
sx q[3];
rz(-2.8521027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.745568) q[0];
sx q[0];
rz(-0.73306274) q[0];
sx q[0];
rz(-2.720604) q[0];
rz(1.1517395) q[1];
sx q[1];
rz(-2.9321509) q[1];
sx q[1];
rz(0.71796012) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67506266) q[0];
sx q[0];
rz(-0.27846876) q[0];
sx q[0];
rz(-1.4896926) q[0];
x q[1];
rz(1.1748208) q[2];
sx q[2];
rz(-0.91063344) q[2];
sx q[2];
rz(-0.3044314) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.3612671) q[1];
sx q[1];
rz(-2.8034489) q[1];
sx q[1];
rz(-1.9664155) q[1];
rz(-pi) q[2];
rz(2.0207383) q[3];
sx q[3];
rz(-1.8186343) q[3];
sx q[3];
rz(-2.7155587) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.026160508) q[2];
sx q[2];
rz(-0.60439503) q[2];
sx q[2];
rz(1.3755414) q[2];
rz(-0.17624217) q[3];
sx q[3];
rz(-0.79564375) q[3];
sx q[3];
rz(-3.0680883) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3714704) q[0];
sx q[0];
rz(-1.4513411) q[0];
sx q[0];
rz(2.0392192) q[0];
rz(0.86028987) q[1];
sx q[1];
rz(-1.5033036) q[1];
sx q[1];
rz(-1.1183429) q[1];
rz(1.5085718) q[2];
sx q[2];
rz(-0.80640651) q[2];
sx q[2];
rz(-2.4666532) q[2];
rz(1.1841487) q[3];
sx q[3];
rz(-1.2332698) q[3];
sx q[3];
rz(-1.5982066) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
