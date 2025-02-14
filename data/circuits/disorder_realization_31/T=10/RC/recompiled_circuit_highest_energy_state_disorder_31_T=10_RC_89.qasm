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
rz(1.0072768) q[0];
sx q[0];
rz(-0.64581031) q[0];
sx q[0];
rz(0.3663775) q[0];
rz(-2.3925048) q[1];
sx q[1];
rz(-1.6269416) q[1];
sx q[1];
rz(0.14263022) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3245904) q[0];
sx q[0];
rz(-2.3329371) q[0];
sx q[0];
rz(-2.3456744) q[0];
rz(-0.5515302) q[2];
sx q[2];
rz(-2.8797134) q[2];
sx q[2];
rz(-0.33360815) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.2097358) q[1];
sx q[1];
rz(-1.0967347) q[1];
sx q[1];
rz(3.0074688) q[1];
x q[2];
rz(-2.5589988) q[3];
sx q[3];
rz(-2.5911463) q[3];
sx q[3];
rz(-0.96874505) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.7819405) q[2];
sx q[2];
rz(-2.6617229) q[2];
sx q[2];
rz(1.5521607) q[2];
rz(1.747067) q[3];
sx q[3];
rz(-0.37140578) q[3];
sx q[3];
rz(2.4594405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8438016) q[0];
sx q[0];
rz(-0.60972917) q[0];
sx q[0];
rz(-0.30526701) q[0];
rz(-1.3097395) q[1];
sx q[1];
rz(-2.7204308) q[1];
sx q[1];
rz(-0.052068204) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.16927707) q[0];
sx q[0];
rz(-2.9575363) q[0];
sx q[0];
rz(0.64223023) q[0];
rz(-2.7676653) q[2];
sx q[2];
rz(-1.2902593) q[2];
sx q[2];
rz(-0.032023059) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.7050208) q[1];
sx q[1];
rz(-1.8364882) q[1];
sx q[1];
rz(2.9712241) q[1];
rz(-0.90272119) q[3];
sx q[3];
rz(-1.4666712) q[3];
sx q[3];
rz(-1.6834761) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.084595844) q[2];
sx q[2];
rz(-1.3330385) q[2];
sx q[2];
rz(-1.4196654) q[2];
rz(2.2325884) q[3];
sx q[3];
rz(-0.82908583) q[3];
sx q[3];
rz(1.7349294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97515714) q[0];
sx q[0];
rz(-1.9571914) q[0];
sx q[0];
rz(-1.3982406) q[0];
rz(-2.1947529) q[1];
sx q[1];
rz(-2.0473174) q[1];
sx q[1];
rz(1.2907226) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9256482) q[0];
sx q[0];
rz(-0.15488347) q[0];
sx q[0];
rz(1.0607222) q[0];
rz(-pi) q[1];
rz(-1.090858) q[2];
sx q[2];
rz(-0.84448123) q[2];
sx q[2];
rz(1.2588071) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-3.0086338) q[1];
sx q[1];
rz(-1.6375638) q[1];
sx q[1];
rz(-3.1125599) q[1];
rz(1.7730784) q[3];
sx q[3];
rz(-1.1014538) q[3];
sx q[3];
rz(0.795006) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.1281841) q[2];
sx q[2];
rz(-0.42542294) q[2];
sx q[2];
rz(-2.181633) q[2];
rz(0.32402447) q[3];
sx q[3];
rz(-0.5158546) q[3];
sx q[3];
rz(0.44017756) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0078761) q[0];
sx q[0];
rz(-2.3634) q[0];
sx q[0];
rz(2.9414951) q[0];
rz(-2.5307185) q[1];
sx q[1];
rz(-0.69823825) q[1];
sx q[1];
rz(-0.79455882) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.824423) q[0];
sx q[0];
rz(-0.18704913) q[0];
sx q[0];
rz(-2.4991799) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7951387) q[2];
sx q[2];
rz(-0.68549978) q[2];
sx q[2];
rz(0.81762527) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.7979203) q[1];
sx q[1];
rz(-1.8842375) q[1];
sx q[1];
rz(-0.15943751) q[1];
x q[2];
rz(-2.4684687) q[3];
sx q[3];
rz(-2.2930721) q[3];
sx q[3];
rz(-2.9150042) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.2574629) q[2];
sx q[2];
rz(-2.5778759) q[2];
sx q[2];
rz(1.5719315) q[2];
rz(-1.3563159) q[3];
sx q[3];
rz(-1.1607728) q[3];
sx q[3];
rz(-3.0089488) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86303478) q[0];
sx q[0];
rz(-0.25535169) q[0];
sx q[0];
rz(0.72364664) q[0];
rz(3.0408995) q[1];
sx q[1];
rz(-2.0932525) q[1];
sx q[1];
rz(0.06631276) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1785051) q[0];
sx q[0];
rz(-1.5076935) q[0];
sx q[0];
rz(-1.0020688) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2264549) q[2];
sx q[2];
rz(-2.0210638) q[2];
sx q[2];
rz(1.7129667) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.7209172) q[1];
sx q[1];
rz(-0.92871341) q[1];
sx q[1];
rz(-0.20963078) q[1];
x q[2];
rz(-2.4661786) q[3];
sx q[3];
rz(-2.6468122) q[3];
sx q[3];
rz(0.29008415) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.940332) q[2];
sx q[2];
rz(-1.2137493) q[2];
sx q[2];
rz(0.21855375) q[2];
rz(2.7471733) q[3];
sx q[3];
rz(-0.68587488) q[3];
sx q[3];
rz(-2.7421537) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0938996) q[0];
sx q[0];
rz(-2.221929) q[0];
sx q[0];
rz(3.1374875) q[0];
rz(0.63402367) q[1];
sx q[1];
rz(-1.5019633) q[1];
sx q[1];
rz(0.9563458) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.078236899) q[0];
sx q[0];
rz(-0.18785297) q[0];
sx q[0];
rz(-2.9211723) q[0];
x q[1];
rz(2.6317503) q[2];
sx q[2];
rz(-1.9467453) q[2];
sx q[2];
rz(-2.816566) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.2507657) q[1];
sx q[1];
rz(-0.27784744) q[1];
sx q[1];
rz(0.020093159) q[1];
rz(1.6878674) q[3];
sx q[3];
rz(-1.9170951) q[3];
sx q[3];
rz(-2.9786144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.52028209) q[2];
sx q[2];
rz(-1.4352398) q[2];
sx q[2];
rz(-2.3069646) q[2];
rz(-2.7708715) q[3];
sx q[3];
rz(-0.70086896) q[3];
sx q[3];
rz(-0.71681517) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
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
rz(-1.9653559) q[0];
sx q[0];
rz(-0.014572425) q[0];
sx q[0];
rz(-0.45005774) q[0];
rz(-2.693148) q[1];
sx q[1];
rz(-2.3095755) q[1];
sx q[1];
rz(-2.4041596) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.099235155) q[0];
sx q[0];
rz(-2.6747534) q[0];
sx q[0];
rz(1.1506266) q[0];
x q[1];
rz(1.7121404) q[2];
sx q[2];
rz(-1.6108043) q[2];
sx q[2];
rz(1.9136946) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.9818283) q[1];
sx q[1];
rz(-2.6609586) q[1];
sx q[1];
rz(0.7483866) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.27100631) q[3];
sx q[3];
rz(-0.79870634) q[3];
sx q[3];
rz(-0.17548156) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.2927148) q[2];
sx q[2];
rz(-2.541031) q[2];
sx q[2];
rz(2.4265477) q[2];
rz(-2.6333366) q[3];
sx q[3];
rz(-1.6786989) q[3];
sx q[3];
rz(1.7432632) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(-pi/2) q[3];
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
rz(2.7364863) q[0];
sx q[0];
rz(-2.5953601) q[0];
sx q[0];
rz(2.7832094) q[0];
rz(-0.37250039) q[1];
sx q[1];
rz(-1.3595711) q[1];
sx q[1];
rz(2.701766) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7067703) q[0];
sx q[0];
rz(-0.79688886) q[0];
sx q[0];
rz(-0.057698542) q[0];
rz(-pi) q[1];
rz(2.7372375) q[2];
sx q[2];
rz(-1.99843) q[2];
sx q[2];
rz(-2.9831487) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.32268128) q[1];
sx q[1];
rz(-2.8214294) q[1];
sx q[1];
rz(2.5500831) q[1];
rz(-pi) q[2];
rz(-1.6251011) q[3];
sx q[3];
rz(-1.2479932) q[3];
sx q[3];
rz(-1.0564547) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0516574) q[2];
sx q[2];
rz(-3.0445485) q[2];
sx q[2];
rz(-0.70866054) q[2];
rz(0.25157252) q[3];
sx q[3];
rz(-2.1422062) q[3];
sx q[3];
rz(-0.034041762) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2947023) q[0];
sx q[0];
rz(-1.0346233) q[0];
sx q[0];
rz(2.5501472) q[0];
rz(-0.017631831) q[1];
sx q[1];
rz(-2.089274) q[1];
sx q[1];
rz(0.10094053) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2786733) q[0];
sx q[0];
rz(-1.8846719) q[0];
sx q[0];
rz(-0.14213965) q[0];
rz(-pi) q[1];
rz(-3.0885126) q[2];
sx q[2];
rz(-0.63737042) q[2];
sx q[2];
rz(-0.87634477) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.7314607) q[1];
sx q[1];
rz(-1.6900151) q[1];
sx q[1];
rz(1.340926) q[1];
x q[2];
rz(-2.9566393) q[3];
sx q[3];
rz(-0.89941632) q[3];
sx q[3];
rz(3.1175478) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.65323222) q[2];
sx q[2];
rz(-2.4182352) q[2];
sx q[2];
rz(1.4867268) q[2];
rz(2.9871353) q[3];
sx q[3];
rz(-1.1052701) q[3];
sx q[3];
rz(2.7753593) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79009295) q[0];
sx q[0];
rz(-0.14227754) q[0];
sx q[0];
rz(-2.067814) q[0];
rz(-1.0723266) q[1];
sx q[1];
rz(-2.8876979) q[1];
sx q[1];
rz(-3.0825739) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8813777) q[0];
sx q[0];
rz(-0.20199595) q[0];
sx q[0];
rz(-2.9750573) q[0];
rz(-1.121465) q[2];
sx q[2];
rz(-2.5856433) q[2];
sx q[2];
rz(-2.1722171) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.7714951) q[1];
sx q[1];
rz(-2.1019249) q[1];
sx q[1];
rz(-2.4350776) q[1];
x q[2];
rz(-2.288184) q[3];
sx q[3];
rz(-1.0685669) q[3];
sx q[3];
rz(0.870992) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.8129639) q[2];
sx q[2];
rz(-1.8495411) q[2];
sx q[2];
rz(0.22895075) q[2];
rz(1.1936584) q[3];
sx q[3];
rz(-2.8173859) q[3];
sx q[3];
rz(-1.4540023) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0070294587) q[0];
sx q[0];
rz(-0.80637359) q[0];
sx q[0];
rz(-0.73643186) q[0];
rz(-1.7560584) q[1];
sx q[1];
rz(-1.7693188) q[1];
sx q[1];
rz(1.7973695) q[1];
rz(-1.6177542) q[2];
sx q[2];
rz(-0.53971077) q[2];
sx q[2];
rz(-2.5765606) q[2];
rz(-1.368461) q[3];
sx q[3];
rz(-2.1740395) q[3];
sx q[3];
rz(-2.5185829) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
