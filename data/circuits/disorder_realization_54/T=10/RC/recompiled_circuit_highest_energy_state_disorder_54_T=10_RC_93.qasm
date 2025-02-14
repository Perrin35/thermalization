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
rz(-0.71159166) q[0];
sx q[0];
rz(1.1049668) q[0];
sx q[0];
rz(8.2850716) q[0];
rz(2.3501514) q[1];
sx q[1];
rz(-2.1005519) q[1];
sx q[1];
rz(-2.0120373) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72973824) q[0];
sx q[0];
rz(-1.751873) q[0];
sx q[0];
rz(-0.84946655) q[0];
rz(0.036058124) q[2];
sx q[2];
rz(-2.4781422) q[2];
sx q[2];
rz(-2.5306338) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.754555) q[1];
sx q[1];
rz(-2.3554152) q[1];
sx q[1];
rz(-1.1494067) q[1];
x q[2];
rz(0.33282354) q[3];
sx q[3];
rz(-1.6037919) q[3];
sx q[3];
rz(0.28256114) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.0396042) q[2];
sx q[2];
rz(-2.3309989) q[2];
sx q[2];
rz(-0.90768901) q[2];
rz(0.11861079) q[3];
sx q[3];
rz(-1.5396996) q[3];
sx q[3];
rz(-2.8906726) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94659656) q[0];
sx q[0];
rz(-2.4793766) q[0];
sx q[0];
rz(-0.10435852) q[0];
rz(1.9508427) q[1];
sx q[1];
rz(-1.933681) q[1];
sx q[1];
rz(-0.45713919) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32923082) q[0];
sx q[0];
rz(-2.9456249) q[0];
sx q[0];
rz(1.4120483) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6847849) q[2];
sx q[2];
rz(-0.94508445) q[2];
sx q[2];
rz(2.5018951) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.10875952) q[1];
sx q[1];
rz(-1.6656414) q[1];
sx q[1];
rz(-1.8873358) q[1];
rz(-pi) q[2];
rz(0.13808226) q[3];
sx q[3];
rz(-0.90681091) q[3];
sx q[3];
rz(-2.9241204) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.55598688) q[2];
sx q[2];
rz(-1.9827236) q[2];
sx q[2];
rz(-2.9846094) q[2];
rz(-2.6202776) q[3];
sx q[3];
rz(-2.3020404) q[3];
sx q[3];
rz(3.1275911) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-0.74612015) q[0];
sx q[0];
rz(-2.5465901) q[0];
sx q[0];
rz(0.34573063) q[0];
rz(-1.455447) q[1];
sx q[1];
rz(-1.2724178) q[1];
sx q[1];
rz(-1.5706496) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.60590505) q[0];
sx q[0];
rz(-2.5890124) q[0];
sx q[0];
rz(1.4009757) q[0];
rz(2.1934671) q[2];
sx q[2];
rz(-1.5678763) q[2];
sx q[2];
rz(2.7013283) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.8040672) q[1];
sx q[1];
rz(-0.21906549) q[1];
sx q[1];
rz(2.9143224) q[1];
x q[2];
rz(0.9216347) q[3];
sx q[3];
rz(-1.4086257) q[3];
sx q[3];
rz(2.0258486) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3674783) q[2];
sx q[2];
rz(-2.9134637) q[2];
sx q[2];
rz(2.1841614) q[2];
rz(1.1227603) q[3];
sx q[3];
rz(-1.9072396) q[3];
sx q[3];
rz(1.3721589) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(0.51266176) q[0];
sx q[0];
rz(-1.4224195) q[0];
sx q[0];
rz(-1.5439532) q[0];
rz(-1.3900025) q[1];
sx q[1];
rz(-1.3163687) q[1];
sx q[1];
rz(-1.66473) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47361237) q[0];
sx q[0];
rz(-1.2897217) q[0];
sx q[0];
rz(-2.5949508) q[0];
rz(-pi) q[1];
rz(1.3324758) q[2];
sx q[2];
rz(-1.4891948) q[2];
sx q[2];
rz(-1.9072744) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.37472607) q[1];
sx q[1];
rz(-1.0488759) q[1];
sx q[1];
rz(0.38066997) q[1];
rz(2.4493205) q[3];
sx q[3];
rz(-0.41818968) q[3];
sx q[3];
rz(-1.632477) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.59820286) q[2];
sx q[2];
rz(-1.2849839) q[2];
sx q[2];
rz(-2.1447935) q[2];
rz(-1.2654842) q[3];
sx q[3];
rz(-0.71804738) q[3];
sx q[3];
rz(-0.27749458) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0811512) q[0];
sx q[0];
rz(-1.5716946) q[0];
sx q[0];
rz(-1.0720217) q[0];
rz(0.95442665) q[1];
sx q[1];
rz(-1.9642893) q[1];
sx q[1];
rz(-1.4929474) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7537751) q[0];
sx q[0];
rz(-0.79083453) q[0];
sx q[0];
rz(-2.8641939) q[0];
rz(-pi) q[1];
rz(0.49333879) q[2];
sx q[2];
rz(-0.78560053) q[2];
sx q[2];
rz(1.0935022) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.49661885) q[1];
sx q[1];
rz(-0.3589093) q[1];
sx q[1];
rz(-0.41618698) q[1];
rz(-pi) q[2];
rz(2.8478283) q[3];
sx q[3];
rz(-0.33393327) q[3];
sx q[3];
rz(1.6373843) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2299049) q[2];
sx q[2];
rz(-0.3868843) q[2];
sx q[2];
rz(1.9913199) q[2];
rz(3.1005499) q[3];
sx q[3];
rz(-1.5464455) q[3];
sx q[3];
rz(2.7400147) q[3];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3732442) q[0];
sx q[0];
rz(-2.0337489) q[0];
sx q[0];
rz(1.9819697) q[0];
rz(1.4609963) q[1];
sx q[1];
rz(-2.9407839) q[1];
sx q[1];
rz(1.9083091) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.999812) q[0];
sx q[0];
rz(-1.8798057) q[0];
sx q[0];
rz(-0.76844333) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.057282863) q[2];
sx q[2];
rz(-1.6241239) q[2];
sx q[2];
rz(-1.48207) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.4602051) q[1];
sx q[1];
rz(-0.80599313) q[1];
sx q[1];
rz(0.56350033) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.308611) q[3];
sx q[3];
rz(-1.0592469) q[3];
sx q[3];
rz(-1.2894323) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.636574) q[2];
sx q[2];
rz(-0.36618149) q[2];
sx q[2];
rz(1.9912857) q[2];
rz(0.73927528) q[3];
sx q[3];
rz(-1.8940247) q[3];
sx q[3];
rz(2.6967743) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66733852) q[0];
sx q[0];
rz(-2.4485454) q[0];
sx q[0];
rz(0.44878238) q[0];
rz(-1.5397286) q[1];
sx q[1];
rz(-2.0346784) q[1];
sx q[1];
rz(1.9805699) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0036573) q[0];
sx q[0];
rz(-1.1627534) q[0];
sx q[0];
rz(-1.28027) q[0];
x q[1];
rz(-2.3660701) q[2];
sx q[2];
rz(-0.65185129) q[2];
sx q[2];
rz(1.3339804) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.4515069) q[1];
sx q[1];
rz(-1.2770318) q[1];
sx q[1];
rz(-0.34505941) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.26412873) q[3];
sx q[3];
rz(-1.0326516) q[3];
sx q[3];
rz(-1.2387741) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.9728969) q[2];
sx q[2];
rz(-1.3316414) q[2];
sx q[2];
rz(1.9608344) q[2];
rz(2.0598038) q[3];
sx q[3];
rz(-1.6094145) q[3];
sx q[3];
rz(0.42657524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5931554) q[0];
sx q[0];
rz(-1.8226382) q[0];
sx q[0];
rz(-2.4539808) q[0];
rz(-1.9718735) q[1];
sx q[1];
rz(-2.1673188) q[1];
sx q[1];
rz(2.7489472) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8075645) q[0];
sx q[0];
rz(-2.0453718) q[0];
sx q[0];
rz(2.9827098) q[0];
rz(-1.7101426) q[2];
sx q[2];
rz(-1.8799439) q[2];
sx q[2];
rz(1.6935371) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.029812977) q[1];
sx q[1];
rz(-1.6583484) q[1];
sx q[1];
rz(-0.47595166) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6527376) q[3];
sx q[3];
rz(-1.382516) q[3];
sx q[3];
rz(-0.57037607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.9795867) q[2];
sx q[2];
rz(-1.0065099) q[2];
sx q[2];
rz(-1.6458192) q[2];
rz(1.5227854) q[3];
sx q[3];
rz(-2.3706172) q[3];
sx q[3];
rz(-2.5405367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
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
rz(2.9252121) q[0];
sx q[0];
rz(-2.7583211) q[0];
sx q[0];
rz(-1.8076757) q[0];
rz(-1.9981617) q[1];
sx q[1];
rz(-0.83088487) q[1];
sx q[1];
rz(0.7935895) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.586602) q[0];
sx q[0];
rz(-1.7287917) q[0];
sx q[0];
rz(-1.642307) q[0];
rz(-pi) q[1];
rz(-0.066448224) q[2];
sx q[2];
rz(-2.1796103) q[2];
sx q[2];
rz(-2.7087536) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.55793437) q[1];
sx q[1];
rz(-2.0730632) q[1];
sx q[1];
rz(1.7947547) q[1];
x q[2];
rz(2.0212444) q[3];
sx q[3];
rz(-1.4739047) q[3];
sx q[3];
rz(0.052122744) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.44635832) q[2];
sx q[2];
rz(-2.0909205) q[2];
sx q[2];
rz(1.0515155) q[2];
rz(2.4255883) q[3];
sx q[3];
rz(-0.99596888) q[3];
sx q[3];
rz(-2.9314465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7150772) q[0];
sx q[0];
rz(-2.2861013) q[0];
sx q[0];
rz(-0.18095428) q[0];
rz(-1.9521693) q[1];
sx q[1];
rz(-1.1782497) q[1];
sx q[1];
rz(-2.1612371) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0553833) q[0];
sx q[0];
rz(-2.6329941) q[0];
sx q[0];
rz(0.17091708) q[0];
x q[1];
rz(-2.0043401) q[2];
sx q[2];
rz(-1.989365) q[2];
sx q[2];
rz(0.43805447) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.7661383) q[1];
sx q[1];
rz(-2.3050637) q[1];
sx q[1];
rz(1.3260654) q[1];
rz(-0.3950903) q[3];
sx q[3];
rz(-1.7655585) q[3];
sx q[3];
rz(2.2421746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.02562) q[2];
sx q[2];
rz(-0.16416922) q[2];
sx q[2];
rz(-1.3244965) q[2];
rz(-2.4888511) q[3];
sx q[3];
rz(-2.1721811) q[3];
sx q[3];
rz(0.038221922) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7272335) q[0];
sx q[0];
rz(-1.4215195) q[0];
sx q[0];
rz(-0.97378578) q[0];
rz(0.7242135) q[1];
sx q[1];
rz(-1.0754633) q[1];
sx q[1];
rz(0.16574688) q[1];
rz(-0.090567055) q[2];
sx q[2];
rz(-1.236785) q[2];
sx q[2];
rz(-1.0192237) q[2];
rz(-2.3694585) q[3];
sx q[3];
rz(-1.6498389) q[3];
sx q[3];
rz(-2.6966358) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
