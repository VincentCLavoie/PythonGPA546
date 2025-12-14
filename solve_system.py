# ==============================================================================
# SOLVE AND DISPLAY RESULTS AT THE TOP (ENHANCED SYSTEM SOLVER)
# ==============================================================================

with solution_box:
    if solve_clicked:
        if not selected_eq_indices:
            st.warning("Select at least one element-wise equation to solve from.")
        elif not selected_vars:
            st.warning("Select at least one variable to solve for in the sidebar.")
        else:
            st.markdown("### Solutions")
            with st.container(border=True):
                # Collect non-trivial selected equations
                selected_eqs = []
                for eq_idx in selected_eq_indices:
                    eq = scalar_eqs[eq_idx]
                    i, j = divmod(eq_idx, 4)
                    if eq is True or eq == True:
                        st.write(f"Equation ({i+1},{j+1}) is identically true; skipping.")
                        continue
                    selected_eqs.append((eq_idx, eq))

                if not selected_eqs:
                    st.write("No non-trivial equations selected.")
                else:
                    # Use all selected equations together as a system
                    eq_list = [eq for _, eq in selected_eqs]
                    
                    # Display which equations are being used
                    eq_indices_str = ", ".join(
                        f"({divmod(eq_idx, 4)[0]+1},{divmod(eq_idx, 4)[1]+1})"
                        for eq_idx, _ in selected_eqs
                    )
                    st.markdown(f"**System of {len(eq_list)} equation(s):** {eq_indices_str}")
                    st.divider()

                    for target in selected_vars:
                        sol = sp.solve(eq_list, target, dict=True)
                        if sol:
                            # If multiple solutions, show all (or first one by default)
                            st.latex(
                                sp.latex(target)
                                + " = "
                                + latex_robotic(sol[0][target])
                            )
                            # Optional: if there are multiple solutions, show count
                            if len(sol) > 1:
                                st.caption(f"(Solution 1 of {len(sol)})")
                        else:
                            st.write(
                                f"❌ No solution for {sp.latex(target)} "
                                f"using system of equations {eq_indices_str}."
                            )
