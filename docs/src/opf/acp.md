# AC-OPF


## Mathematical Formulation

The ACOPF model considered in PGLearn is presented below.

```math
\begin{align}
    \min_{\PG, \QG, \PF, \QF, \PT, \QT, \VM, \VA} \quad
    & \sum_{i \in \NODES} \sum_{j \in \GENERATORS_{i}} c_j \PG_j + c_0 \label{model:acopf:obj} \\
    \text{s.t.} \quad \quad \quad
    & \sum_{j\in\GENERATORS_i}\PG_j - \sum_{j\in\LOADS_i}\PD_j - \GS_i \VM_i^2 = \sum_{e \in \mathcal{E}^{+}_{i}}  \PF_{e} + \sum_{e \in \mathcal{E}^{-}_{i}} \PT_{e}
        & \forall i &\in \NODES
        \label{model:acopf:kirchhoff:active} \\
    & \sum_{j\in\GENERATORS_i}\QG_j - \sum_{j\in\LOADS_i}\QD_j + \BS_i \VM_i^2 = \sum_{e \in \mathcal{E}^{+}_{i}} \QF_{e} + \sum_{e \in \mathcal{E}^{-}_{i}} \QT_{e}
        & \forall i &\in \NODES
        \label{model:acopf:kirchhoff:reactive} \\
    & \PF_{e} = \gff_{e}\VM_i^2 + \gft_{e} \VM_i \VM_j \cos(\VA_i-\VA_j) + \bft_{e} \VM_i \VM_j \sin(\VA_i-\VA_j)
        & \forall e = (i,j) &\in \EDGES
        \label{model:acopf:ohm:active:from} \\
    & \QF_{e} = -\bff_{e} \VM_i^2 - \bft_{e}\VM_i \VM_j \cos(\VA_i-\VA_j) + \gft_{e} \VM_i \VM_j \sin(\VA_i-\VA_j)
        & \forall e = (i,j) &\in \EDGES
        \label{model:acopf:ohm:reactive:from} \\
    & \PT_{e} = \gtt_{e}\VM_j^2 + \gtf_{e} \VM_i \VM_j \cos(\VA_i-\VA_j) - \btf_{e} \VM_i \VM_j \sin(\VA_i-\VA_j)
        & \forall e = (i,j) &\in \EDGES
        \label{model:acopf:ohm:active:to} \\
    & \QT_{e} = -\btt_{e} \VM_j^2 - \btf_{e}\VM_i \VM_j \cos(\VA_i-\VA_j) - \gtf_{e} \VM_i \VM_j \sin(\VA_i-\VA_j)
        & \forall e = (i,j) &\in \EDGES
        \label{model:acopf:ohm:reactive:to} \\
    & (\PF_{e})^2 + (\QF_{e})^2 \leq \overline{S}_{e}^2
        & \forall e &\in \EDGES
        \label{model:acopf:thrmbound:from} \\
    & (\PT_{e})^2 + (\QT_{e})^2 \leq \overline{S}_{e}^2
        & \forall e &\in \EDGES
        \label{model:acopf:thrmbound:to} \\
    & \dvamin_{e} \leq \VA_i - \VA_j \leq \dvamax_{e}
        & \forall e = (i,j) &\in \EDGES
        \label{model:acopf:angledifference} \\
    & \VA_\text{ref} = 0
        \label{model:acopf:slackbus} \\
    & \pgmin_{i} \leq \PG_i \leq \pgmax_{i}
        & \forall i &\in \GENERATORS
        \label{model:acopf:pgbound} \\
    & \qgmin_{i} \leq \QG_i \leq \qgmax_{i}
        & \forall i &\in \GENERATORS
        \label{model:acopf:qgbound} \\
    & \vmmin_{i} \leq \VM_i \leq \vmmax_{i}
        & \forall i &\in \NODES
        \label{model:acopf:vmbound} \\
    & {-}\overline{S}_{e} \leq  \PF_{e} \leq \overline{S}_{e}
        & \forall e &\in \EDGES
        \label{model:acopf:pfbound} \\
    & {-}\overline{S}_{e} \leq  \QF_{e} \leq \overline{S}_{e}
        & \forall e &\in \EDGES
        \label{model:acopf:qfbound} \\
    & {-}\overline{S}_{e} \leq  \PT_{e} \leq \overline{S}_{e}
        & \forall e &\in \EDGES
        \label{model:acopf:ptbound} \\
    & {-}\overline{S}_{e} \leq  \QT_{e} \leq \overline{S}_{e}
        & \forall e &\in \EDGES
        \label{model:acopf:qtbound}
\end{align}
```


### Variables

* ``\VM \in \mathbb{R}^{N}``: nodal voltage magnitude
* ``\VA \in \mathbb{R}^{N}``: nodal voltage angle
* ``\PG \in \mathbb{R}^{G}``: active power dispatch
* ``\QG \in \mathbb{R}^{G}``: reactive power dispatch
* ``\PF \in \mathbb{R}^{E}``: active power flow "from"
* ``\QF \in \mathbb{R}^{E}``: reactive power flow "from"
* ``\PT \in \mathbb{R}^{E}``: active power flow "to"
* ``\QT \in \mathbb{R}^{E}``: reactive power flow "to"

### Objective

The objective function ``\eqref{model:acopf:obj}`` minimizes the cost of active power generation.

!!! todo
    PGLearn currently supports only linear cost functions.
    Support for quadratic functions is planned for a later stage; please open an issue if 
    you would like to request this feature.

### Constraints

* ``\eqref{model:acopf:kirchhoff:active}-\eqref{model:acopf:kirchhoff:reactive}``:
    active and reactive Kirchhoff's current law.
* ``\eqref{model:acopf:ohm:active:from}-\eqref{model:acopf:ohm:reactive:to}``:
    Ohm's law expressing power flows as a function of nodal voltages.
* ``\eqref{model:acopf:thrmbound:from}-\eqref{model:acopf:thrmbound:to}``: Thermal limits
* ``\eqref{model:acopf:angledifference}``: Voltage angle deviation constraints.
* ``\eqref{model:acopf:slackbus}``: this constraint fixes the voltage angle of the 
    reference (slack) bus to zero.
* ``\eqref{model:acopf:pgbound}-\eqref{model:acopf:qgbound}``: Active/reactive generation limits
* ``\eqref{model:acopf:vmbound}``: Nodal voltage magnitude limits
* ``\eqref{model:acopf:pfbound}-\eqref{model:acopf:qtbound}``: Power flow bounds, derived from thermal limits

!!! info
    Although power flow variable bounds ``\eqref{model:acopf:pfbound}-\eqref{model:acopf:qtbound}``
    are redundant with thermal limits ``\eqref{model:acopf:thrmbound:from}-\eqref{model:acopf:thrmbound:to}``, 
    their inclusion improves the performance of interior-point solvers like Ipopt.


## Data format

### Primal solution

| Variable  | Data | Size  | Description 
|:----------|:-----|:------|:----------------------------------|
|  ``\VM``  | `vm` | ``N`` | Nodal voltage magnitude
|  ``\VA``  | `va` | ``N`` | Nodal voltage angle
|  ``\PG``  | `pg` | ``G`` | Active power generation
|  ``\QG``  | `pg` | ``G`` | Reactive power generation
|  ``\PF``  | `pf` | ``E`` | Active power flow (from)
|  ``\PT``  | `pt` | ``E`` | Active power flow (to)
|  ``\QF``  | `qf` | ``E`` | Reactive power flow (from)
|  ``\QT``  | `qt` | ``E`` | Reactive power flow (to)

### Dual solution

| Constraint                                   | Data         | Size  |
|:---------------------------------------------|:-------------|:------|
| ``\eqref{model:acopf:kirchhoff:active}``     | `kcl_p`      | ``N`` |
| ``\eqref{model:acopf:kirchhoff:reactive}``   | `kcl_q`      | ``N`` |
| ``\eqref{model:acopf:ohm:active:from}``      | `ohm_pf`     | ``E`` |
| ``\eqref{model:acopf:ohm:reactive:from}``    | `ohm_qf`     | ``E`` |
| ``\eqref{model:acopf:ohm:active:to}``        | `ohm_pt`     | ``E`` |
| ``\eqref{model:acopf:ohm:reactive:to}``      | `ohm_qt`     | ``E`` |
| ``\eqref{model:acopf:thrmbound:from}``       | `sm_fr`      | ``E`` |
| ``\eqref{model:acopf:thrmbound:to}``         | `sm_to`      | ``E`` |
| ``\eqref{model:acopf:angledifference}``      | `va_diff`    | ``E`` |
| ``\eqref{model:acopf:slackbus}``             | `slack_bus`  | ``N`` |
| ``\eqref{model:acopf:pgbound}``              | `pg`         | ``G`` |
| ``\eqref{model:acopf:qgbound}``              | `qg`         | ``G`` |
| ``\eqref{model:acopf:vmbound}``              | `vm`         | ``N`` |
| ``\eqref{model:acopf:pfbound}``              | `pf`         | ``E`` |
| ``\eqref{model:acopf:qfbound}``              | `qf`         | ``E`` |
| ``\eqref{model:acopf:ptbound}``              | `pt`         | ``E`` |
| ``\eqref{model:acopf:qtbound}``              | `qt`         | ``E`` |
